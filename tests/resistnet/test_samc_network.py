import copy
import os

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import resistnet
from resistnet.samc_network import ResistanceNetworkSAMC
from resistnet.model_optimisation import ModelRunnerTPE

VARIABLES = ["run_mm_cyr", "tmp_dc_cyr", "dor_pc_pva", "sgr_dk_rav"]

# 60-digit references (dense Gaussian elimination in decimal arithmetic; see
# test_CFPT._decimal_cond_passage) for the conditional passage times from
# each site to one site, with movement weights between reaches from a
# directional variable under the original (unbounded) resistance-to-
# conductance mapping. Plain sparse LU solves were off by up to 230%
# (runoff, to site 9) and by a factor of 1e31 (temperature, to site 5)
REF_RUN_DIR_TO_SITE9 = [
    588.3286525909401, 561.3091226587889, 574.4379504133846,
    421.97731350087105, 147.18458100692195, 145.93435543276166,
    141.13975585634714, 106.32495420995706, 71.70723308275628, 0.0,
    24.697134857097442, 25343.82141073529, 25787.654863361222,
    25826.699635403125
]
REF_TMP_DIR_TO_SITE5 = [
    4.323455642275723e+16, 4.323455642275683e+16, 665.7155102214871,
    63.71551022151208, 53.75531121651189, 0.0, 64.48942420989053,
    4.323455642275689e+16, 4.3234556422756904e+16, 1.6225927682921414e+31,
    1.6225927682921414e+31, 1.6225927682921414e+31, 1.6225927682921432e+31,
    1.6225927682921468e+31
]


def model(spec):
    """spec: {variable: (weight, transformation, shape)} (not directional)"""
    ind = []
    for v in VARIABLES:
        ind += [1, *spec[v], 0] if v in spec else [0, 1.0, 0, 1, 0]
    return ind


@pytest.fixture(scope="module")
def samc_net(tmp_path_factory):
    out = tmp_path_factory.mktemp("samc")
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    return ResistanceNetworkSAMC(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"),
        sizes=os.path.join(data, "nmt.popSizes.txt"),
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=VARIABLES,
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={v: "ARITH" for v in VARIABLES},
        out=str(out / "samc"), verbose=False
    )


def test_outlet_reaches_linked(samc_net, tmp_path):
    # In nmt.gpkg, the network's lowest reach drains into a reach of the
    # layer, which becomes the single root reach
    n_comp, _ = sp.csgraph.connected_components(samc_net._adj,
                                                directed=False)
    assert n_comp == 1
    root = [d[samc_net.reachid_col] for _, _, d in
            samc_net._K.in_edges(samc_net._origin, data=True)]
    assert root == [70648532]

    # Clipped at that reach, two reaches meet at the outlet: they are linked
    # directly, otherwise the network splits in two
    import geopandas as gpd
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    reaches = gpd.read_file(os.path.join(data, "nmt.gpkg"))
    below = dict(zip(reaches["HYRIV_ID"], reaches["NEXT_DOWN"]))
    drop = [70648532]
    while below[drop[-1]] in below:
        drop.append(below[drop[-1]])
    clipped = str(tmp_path / "clipped.gpkg")
    reaches[~reaches["HYRIV_ID"].isin(drop)].to_file(clipped, driver="GPKG")
    net = ResistanceNetworkSAMC(
        network=None, shapefile=clipped, sizes=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=VARIABLES[:1],
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={VARIABLES[0]: "ARITH"}, out=str(tmp_path / "clipped"),
        verbose=False)
    n_comp, _ = sp.csgraph.connected_components(net._adj, directed=False)
    assert n_comp == 1
    assert len(net._edge_site_indices) == 14
    assert len(net._site_states) == 14
    idx = {e: i for i, e in enumerate(net._edge_order)}
    root = sorted(idx[d[net.reachid_col]] for _, _, d in
                  net._K.in_edges(net._origin, data=True))
    assert sorted(net._edge_order[i] for i in root) == [70648307, 70648308]
    assert net._adj[root[0], root[1]] == -1
    assert net._adj[root[1], root[0]] == -1
    # On the lattice, the two reaches meet at the outlet node
    lat = net._lattice
    outlet = lat.states([net._origin])[0]
    assert lat.pos[outlet] == 0
    assert (lat.b[lat.a == outlet] >= len(lat.nodes)).sum() == 2


def test_lattice_sites_and_absorption(samc_net):
    # Each site is the state of its snapped node, absorbed with
    # probability 1/(2 Ne) per visit; nothing else absorbs
    lat = samc_net._lattice
    sites = list(samc_net._points_snapped.keys())
    assert samc_net._site_states == lat.states(sites)
    R = samc_net._R_lat
    np.testing.assert_allclose(
        R[samc_net._site_states],
        [1 / (2 * samc_net._sizes[s]) for s in sites])
    assert np.count_nonzero(R) == len(sites)
    # About one state per km of river (plus the nodes)
    km = sum(d[samc_net.length_col] for _, _, d in
             samc_net._K.edges(data=True))
    np.testing.assert_allclose(lat.length.sum() / 2, km)
    assert abs(lat.n - len(lat.nodes) - (km - len(samc_net._edge_order))) \
        < 0.5 * len(samc_net._edge_order)


@pytest.mark.parametrize("var,site,ref", [
    ("run_mm_cyr", 9, REF_RUN_DIR_TO_SITE9),
    ("tmp_dc_cyr", 5, REF_TMP_DIR_TO_SITE5)])
def test_cfpt_precision(samc_net, var, site, ref):
    # Passage times on extremely uneven moves between reaches (a numerical
    # test of CFPT)
    from resistnet.CFPT import CFPT, transition_matrix
    from resistnet.utils import minmax, minmax_nonzero
    v = samc_net._predictors[var].to_numpy(dtype=float)
    v = (v - v.min()) / (v.max() - v.min())
    rows, cols = samc_net._adj.nonzero()
    w = minmax_nonzero(1 / minmax_nonzero(minmax(v[cols] - v[rows])))
    W = sp.csr_matrix((w, (rows, cols)), shape=samc_net._adj.shape)
    T = CFPT(transition_matrix(W, samc_net._R), samc_net._edge_site_indices)
    np.testing.assert_allclose(T[:, site], ref, rtol=1e-10)


def test_resistance_is_the_standard_models(samc_net):
    # Each reach's resistance per km is the standard model's (ResistanceGA's
    # surface): 1 + (K - 1) x the weighted mean of the transformed
    # variables; the walker's weights are samc's conductance on the
    # lattice, 2 / ((r_i + r_j) x link length)
    from resistnet.resistance_network import ResistanceNetwork
    ind = model({"run_mm_cyr": (0.7, 3, 10), "sgr_dk_rav": (0.4, 5, 2)})
    for k in [1.0, 37.0, 500.0]:
        rho = samc_net.resistance_surface(ind + [k]).to_numpy()
        std = ResistanceNetwork._build_composite_surface(
            samc_net, ind + [k]).to_numpy()
        np.testing.assert_allclose(rho, 1 + (k - 1) * std, rtol=1e-12)
        assert rho.min() >= 1 and rho.max() <= k + 1e-9
        lat = samc_net._lattice
        r = lat.state_resistance(rho)
        W = samc_net._movement_weights(ind + [k]).tocsr()
        np.testing.assert_allclose(
            np.asarray(W[lat.a, lat.b]).ravel(),
            2 / ((r[lat.a] + r[lat.b]) * lat.length), rtol=1e-12)
    # K = 1 is distance alone: conductance 1 / link length
    W = samc_net._movement_weights(ind + [1.0]).tocsr()
    lat = samc_net._lattice
    np.testing.assert_allclose(np.asarray(W[lat.a, lat.b]).ravel(),
                               1 / lat.length, rtol=1e-12)
    # Without a maximum resistance in the model, the default is used
    assert samc_net.max_resistance(ind) == 100.0


def test_directional_variables_refused(samc_net):
    ind = model({"run_mm_cyr": (1.0, 0, 1)})
    ind[4] = 1
    with pytest.raises(ValueError, match="directional"):
        samc_net._movement_weights(ind + [10.0, 0.0])
    assert samc_net.effective_directions(ind) == {}


def test_evaluate_and_output(samc_net):
    ind = model({"run_mm_cyr": (0.7, 3, 10), "sgr_dk_rav": (0.4, 0, 1)})
    fitness, res = samc_net.evaluate(ind + [30.0, 0.01])
    assert np.isfinite(fitness)
    # loglik, r2m, aic, delta_aic_null (as for the base model)
    assert np.all(np.isfinite(res[:4]))
    assert fitness == res[2]

    cfpt, edges = samc_net.model_output(ind + [30.0, 0.01])
    n = len(samc_net._site_states)
    assert cfpt.shape == (n, n)
    assert np.all(np.diag(cfpt) == 0)
    assert np.all(cfpt[~np.eye(n, dtype=bool)] >= 1)
    # Each reach's resistance per km, and its middle state's moves
    # downstream and upstream (no absorption there: they sum to 1)
    np.testing.assert_allclose(
        edges[:, 0], samc_net.resistance_surface(ind + [30.0]).to_numpy())
    has = ~np.isnan(edges[:, 1])
    assert has.mean() > 0.5
    np.testing.assert_allclose(edges[has, 1] + edges[has, 2], 1)
    # A positive flow rate makes moving upstream harder
    assert np.all(edges[has, 1] > 0.5)


def test_search_space(samc_net):
    # Variables are never directional on the lattice; the flow term is
    runner = ModelRunnerTPE(samc_net, seed=1, verbose=False)
    runner.set_tpe_parameters(max_shape=10)
    runner.initialize_space()
    from hyperopt.pyll.stochastic import sample
    rng = np.random.default_rng(1)
    draws = [sample(runner.space, rng=rng) for _ in range(300)]
    assert {d["run_mm_cyr_asym"] for d in draws} == {0}
    assert {d["flow_dir"] for d in draws} == {-1, 0, 1}
    rates = [d["flow_rate"] for d in draws]
    assert min(rates) >= 1e-3 and max(rates) <= samc_net.max_flow
    runner.set_tpe_parameters(max_shape=10, fixFlow=True)
    runner.initialize_space()
    draws = [sample(runner.space, rng=rng) for _ in range(50)]
    assert {d["flow_dir"] for d in draws} == {0}


def test_surface_parameters(samc_net):
    ind = model({"run_mm_cyr": (0.7, 3, 10), "sgr_dk_rav": (0.4, 0, 1)})
    # 2 x (transformation, shape) + 1 relative weight
    assert samc_net.surface_parameters(ind) == 5
    # + the maximum resistance, when it is part of the model
    assert samc_net.surface_parameters(ind + [50.0]) == 6
    net = copy.copy(samc_net)
    net.fix_resist = 50.0
    assert net.surface_parameters(ind + [50.0]) == 5
    # + the flow rate
    assert samc_net.surface_parameters(ind + [50.0, 0.0]) == 6
    assert samc_net.surface_parameters(ind + [50.0, -0.1]) == 7
    # Flow alone: only the rate (no maximum resistance without variables)
    assert samc_net.surface_parameters(model({}) + [50.0, 0.1]) == 1


def test_run_tpe_samc(samc_net, tmp_path):
    runner = ModelRunnerTPE(samc_net, seed=1, verbose=False)
    out = str(tmp_path / "run")
    runner.run_tpe(max_evals=3, threads=1, reps=2, max_shape=10,
                   n_startup=2, out=out, verbose=False, plot=False)
    hof = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t", na_values="-")
    assert all(f"{v}_asym" in hof.columns for v in VARIABLES)
    for v in VARIABLES:
        assert hof[f"{v}_asym"].dropna().isin([0]).all()
    assert hof["max_resistance"].between(1, samc_net.max_resist).all()
    assert hof["flow"].abs().le(samc_net.max_flow).all()
    # AIC differences from the reference models (positive = better)
    null = pd.read_csv(f"{out}.Null-Model.tsv", sep="\t", index_col=0)
    for col, ref in [("delta_aic_distance", "samc_distance"),
                     ("delta_aic_flow", "samc_flow")]:
        np.testing.assert_allclose(
            hof[col], null.loc[ref, "aic"] - hof["aic"], rtol=1e-9)
    assert os.path.exists(f"{out}.Model-Average.CFPTMatrix.tsv")
    edges = pd.read_csv(f"{out}.Model-Average.TransitionEdges.tsv",
                        sep="\t")
    assert list(edges.columns) == ["EDGE_ID", "Resistance", "Downstream",
                                   "Upstream"]
    import json
    with open(f"{out}.modelSpec.json") as fh:
        spec = json.load(fh)
    assert [m["rank"] for m in spec["models"]] == list(
        range(1, len(hof) + 1))
    np.testing.assert_allclose([m["aic"] for m in spec["models"]],
                               hof["aic"])
    assert all(m["cell_km"] == 1.0 for m in spec["models"])


@pytest.mark.parametrize("stored", [("o", "d"), ("d", "o")])
def test_extend_confluence_orientation(stored):
    # The new origin must be the end of the downstream reach away from the
    # current origin, however the edge's endpoints are stored
    import networkx as nx
    net = ResistanceNetworkSAMC.__new__(ResistanceNetworkSAMC)
    net.reachid_col, net.infer_origin = "HYRIV_ID", "NEXT_DOWN"
    net._K = nx.Graph()
    net._K.add_edge("a", "o", HYRIV_ID=1, NEXT_DOWN=3, EDGE_ID=0)
    net._K.add_edge("b", "o", HYRIV_ID=2, NEXT_DOWN=3, EDGE_ID=1)
    net._G = nx.Graph(net._K)
    net._G.add_edge(*stored, HYRIV_ID=3, NEXT_DOWN=0)
    net._origin = "o"
    net.extend_confluence(1, 2)
    assert net._origin == "d"
    assert net._reach_origin == 3
    assert net._K.edges["o", "d"]["EDGE_ID"] == 2


def test_uniform_ne(tmp_path):
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    net = ResistanceNetworkSAMC(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"), sizes=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=VARIABLES[:1], inmat=None,
        agg_opts={VARIABLES[0]: "ARITH"}, out=str(tmp_path / "ne"),
        verbose=False, ne=50, cell_km=2.0
    )
    R = net._R_lat
    np.testing.assert_allclose(R[net._site_states], 1 / (2 * 50))
    assert np.count_nonzero(R) == len(net._site_states)
    # Coarser lattice: links of about 2 km
    assert net._lattice.cell_km == 2.0
    assert np.median(net._lattice.length) > 1.5


def test_null_models(samc_net):
    table = samc_net.evaluate_null_model()
    assert list(table.index) == [
        "distance_only", "samc_distance", "samc_flow", "null"]
    # Distance alone = any model with a maximum resistance of 1
    ind = model({"run_mm_cyr": (1.0, 4, 20)}) + [1.0]
    _, res = samc_net.evaluate(ind)
    assert table.loc["samc_distance", "loglik"] == pytest.approx(res[0])
    # Only the MLPE parameters (intercept, slope, 2 variances) are counted
    assert table.loc["samc_distance", "aic"] == pytest.approx(
        2 * 4 - 2 * res[0])
    # Flow alone: the fitted rate, plus one parameter; it includes distance
    # alone (rate 0), so it fits at least as well
    rate = table.loc["samc_flow", "flow"]
    _, res = samc_net.evaluate(model({}) + [1.0, rate])
    assert table.loc["samc_flow", "loglik"] == pytest.approx(res[0])
    assert table.loc["samc_flow", "aic"] == pytest.approx(2 * 5 - 2 * res[0])
    assert (table.loc["samc_flow", "loglik"]
            >= table.loc["samc_distance", "loglik"])
    for r in [rate * 0.8, rate * 1.25]:
        _, res = samc_net.evaluate(model({}) + [1.0, r])
        assert res[0] <= table.loc["samc_flow", "loglik"] + 1e-6


def test_upstream_changes(samc_net):
    table = samc_net.upstream_changes()
    assert list(table.index) == VARIABLES
    shares = table[["rises", "no_change", "falls"]].astype(float)
    np.testing.assert_allclose(shares.sum(axis=1), 1)
    # Each upstream step goes from a reach's downstream reach to the reach
    # (NEXT_DOWN in the reach layer); the root reach has no downstream
    # reach in the network
    import geopandas as gpd
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    reaches = gpd.read_file(os.path.join(data, "nmt.gpkg")).set_index(
        "HYRIV_ID").loc[list(samc_net._edge_order)]
    linked = reaches["NEXT_DOWN"].isin(reaches.index)
    assert (table["steps"] == linked.sum()).all()
    assert linked.sum() == len(samc_net._edge_order) - 1
    for v in VARIABLES:
        change = (reaches.loc[linked, v].to_numpy()
                  - reaches.loc[reaches.loc[linked, "NEXT_DOWN"], v]
                  .to_numpy())
        np.testing.assert_allclose(
            shares.loc[v], [np.mean(change > 0), np.mean(change == 0),
                            np.mean(change < 0)])
    prefix = samc_net.output_prefix
    assert os.path.exists(f"{prefix}.upstreamChanges.pdf")
    written = pd.read_csv(f"{prefix}.upstreamChanges.tsv", sep="\t",
                          index_col=0)
    np.testing.assert_allclose(written["rises"], shares["rises"])


def _raw(net):
    """Raw values of the predictors, from their 0-10 values and ranges."""
    raw = net._predictors.copy()
    for v, (lo, hi) in net._predictor_range.items():
        raw[v] = lo + raw[v] / 10 * (hi - lo)
    return raw


@pytest.mark.parametrize("t,s", [(0, 1), (1, 3), (2, 5), (3, 4), (4, 6),
                                 (5, 2), (6, 8), (7, 3), (8, 2)])
def test_predict_reproduces_model(samc_net, t, s):
    # Applying the saved specification to the training reaches' raw values
    # gives exactly the fitted model's resistances per km, and the moves
    # from each reach's middle state
    import resistnet.predict as pr
    ind = model({"sgr_dk_rav": (0.6, t, s),
                 "run_mm_cyr": (0.4, 3, 7)}) + [37.0, -0.04]
    spec = samc_net.model_spec(ind)
    assert spec["type"] == "samc" and spec["cell_km"] == 1.0
    reaches = _raw(samc_net)
    reaches.insert(0, "HYRIV_ID", reaches.index.to_numpy())
    length = {d[samc_net.reachid_col]: d[samc_net.length_col]
              for _, _, d in samc_net._K.edges(data=True)}
    reaches["LENGTH_KM"] = [length[e] for e in samc_net._edge_order]
    table, clamped = pr.predict_reaches(
        reaches, {"reachid_col": "HYRIV_ID", "length_col": "LENGTH_KM",
                  "models": [spec]})
    assert all(v == 0 for v in clamped.values())
    np.testing.assert_allclose(
        table["resistance"],
        samc_net.resistance_surface(ind).to_numpy(), rtol=1e-12)
    _, edges = samc_net.model_output(ind)
    np.testing.assert_allclose(table["p_downstream"], edges[:, 1],
                               rtol=1e-12)


def test_predict_script(samc_net, tmp_path):
    # The command-line script writes the predictions for a reach layer
    import geopandas as gpd
    import json
    import subprocess
    import sys
    import resistnet.predict as pr
    ind = model({"sgr_dk_rav": (0.6, 7, 3),
                 "run_mm_cyr": (0.4, 1, 4)}) + [20.0, 0.03]
    spec = {"reachid_col": "HYRIV_ID", "length_col": "LENGTH_KM",
            "models": [samc_net.model_spec(ind)]}
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    reaches = gpd.read_file(os.path.join(data, "nmt.gpkg"))
    reaches = reaches[reaches["HYRIV_ID"].isin(samc_net._edge_order)]
    layer = str(tmp_path / "reaches.gpkg")
    reaches.to_file(layer, driver="GPKG")
    table, _ = pr.predict_reaches(reaches, spec)
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(spec))
    script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                          "predictResistnet.py")
    out = str(tmp_path / "pred")
    p = subprocess.run([sys.executable, script, "-s", layer, "-m",
                        str(path), "-o", out], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    written = pd.read_csv(f"{out}.predicted.tsv", sep="\t")
    np.testing.assert_allclose(written["p_downstream"],
                               table["p_downstream"], rtol=1e-12)
    np.testing.assert_allclose(written["resistance"], table["resistance"],
                               rtol=1e-12)
    assert os.path.exists(f"{out}.predicted.gpkg")


def test_predict_clamps_to_training_range(samc_net):
    import resistnet.predict as pr
    spec = samc_net.model_spec(
        model({"sgr_dk_rav": (1.0, 5, 3)}) + [10.0, 0.0])
    var = spec["variables"]["sgr_dk_rav"]
    lo, hi = var["min"], var["max"]
    inside = pr.variable_values([lo, hi], var)
    outside = pr.variable_values([lo - 5, hi + 5], var)
    np.testing.assert_allclose(outside, inside)
    assert pr.out_of_range([lo - 1, lo, hi, hi + 1], var) == 0.5


def test_flow_term(samc_net):
    lat = samc_net._lattice
    for rate in [0.2, -0.05]:
        # Flow alone: weights are exp(-rate * km moved upstream) / length
        W = samc_net._movement_weights(model({}) + [10.0, rate]).tocsr()
        np.testing.assert_allclose(
            np.asarray(W[lat.a, lat.b]).ravel(),
            np.exp(-rate * lat.upstream) / lat.length, rtol=1e-12)
        # With variables, it multiplies their weights; it is per km, so
        # not affected by the maximum resistance
        spec = model({"run_mm_cyr": (1.0, 0, 1)})
        W0 = samc_net._movement_weights(spec + [10.0, 0.0]).tocsr()
        Wf = samc_net._movement_weights(spec + [10.0, rate]).tocsr()
        np.testing.assert_allclose(
            np.asarray(Wf[lat.a, lat.b]).ravel(),
            np.asarray(W0[lat.a, lat.b]).ravel()
            * np.exp(-rate * lat.upstream), rtol=1e-12)
    # Neither variables nor flow: no model
    assert samc_net._movement_weights(model({}) + [10.0, 0.0]) is None
    fitness, _ = samc_net.evaluate(model({}) + [10.0, 0.0])
    assert fitness == float("-inf")


REFERENCE = os.path.join(os.path.dirname(__file__), "reference")
# Models with references from samc::cond_passage (reference/README.md):
# {variable: (weight, transformation, shape)}, maximum resistance per km
# and flow rate; None = distance alone (resistance 1 per km)
R_SAMC_MODELS = {
    "equal": None,
    "sgr_sym": ({"sgr_dk_rav": (1.0, 0, 1.0)}, 50.0, 0.0),
    "sgr_flow": ({"sgr_dk_rav": (1.0, 5, 3.0)}, 20.0, 0.02),
}


def _ref_net(out):
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    v = ["sgr_dk_rav", "ele_mt_cav"]
    return ResistanceNetworkSAMC(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"),
        sizes=os.path.join(data, "nmt.popSizes.txt"),
        coords=os.path.join(data, "nmt.pointCoords.txt"), variables=v,
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={x: "ARITH" for x in v}, out=os.path.join(str(out), "ref"),
        verbose=False)


@pytest.fixture(scope="module")
def ref_net(tmp_path_factory):
    # The variables of the R reference models
    return _ref_net(tmp_path_factory.mktemp("ref"))


def _ref_transition_matrix(net, name):
    from resistnet.CFPT import transition_matrix
    if R_SAMC_MODELS[name] is None:
        return transition_matrix(net._lattice.weights(), net._R_lat)
    spec, K, rate = R_SAMC_MODELS[name]
    ind = []
    for v in net._predictors.columns:
        ind += [1, *spec[v], 0] if v in spec else [0, 1.0, 0, 1.0, 0]
    return net._transition_matrix(ind + [K, rate])


@pytest.mark.parametrize("name", sorted(R_SAMC_MODELS))
def test_cfpt_matches_r_samc_nmt(ref_net, name):
    # Conditional passage times among the NMT sites on the lattice, as
    # samc::cond_passage (reference/nmt_samc_cond_passage.R)
    from resistnet.CFPT import CFPT
    ref = np.loadtxt(os.path.join(REFERENCE,
                                  f"nmt_samc_lattice_{name}.tsv"))
    Q = _ref_transition_matrix(ref_net, name)
    T = CFPT(Q, ref_net._site_states, ref_net._R_lat)
    off = ~np.eye(len(T), dtype=bool)
    np.testing.assert_allclose(T[off], ref[off], rtol=1e-10)


def test_mlpe_reference_inputs_nmt(ref_net):
    # The genetic distances and the SAMC model distances behind the lme4
    # references (test_MLPE.py) are what ResistNet computes for NMT
    from resistnet.CFPT import CFPT
    Y = np.loadtxt(os.path.join(REFERENCE, "nmt_mlpe_Y.tsv"))
    # (input parsing can differ in the last digit between pandas versions)
    np.testing.assert_allclose(ref_net._gendist, Y, rtol=1e-12, atol=0)
    X = np.loadtxt(os.path.join(REFERENCE, "nmt_mlpe_X_samc_sgr_flow.tsv"))
    T = CFPT(_ref_transition_matrix(ref_net, "sgr_flow"),
             ref_net._site_states, ref_net._R_lat)
    np.testing.assert_allclose((T + T.T) / 2, X, rtol=1e-12)


def _simulate_samc(tmp_path, spec, flow, seed=4, samples=12):
    from resistnet.samc_network import SimResistanceNetworkSAMC
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    spec_file = tmp_path / "spec.tsv"
    spec.to_csv(spec_file, sep="\t", index=False)
    sim = SimResistanceNetworkSAMC(
        network=os.path.join(data, "nmt.network"), reachid_col="HYRIV_ID",
        length_col="LENGTH_KM", seed=seed, max_resistance=20.0, flow=flow)
    out = str(tmp_path / "sim")
    sim.simulate(str(spec_file), 1, samples, out)
    return data, f"{out}_1"


def test_simulate_samc_recovers_flow(tmp_path):
    # Distances simulated under flow alone are recovered by fitting flow
    # alone to them, with the same rate
    spec = pd.DataFrame(columns=["VAR", "WEIGHT", "TRANSFORM", "SHAPE"])
    data, base = _simulate_samc(tmp_path, spec, flow=0.05)
    sim = pd.read_csv(f"{base}.ResistanceMatrix.tsv", sep="\t", index_col=0)
    assert np.allclose(sim.to_numpy(), sim.to_numpy().T)
    assert sim.to_numpy().max() == pytest.approx(1.0)
    cfpt = pd.read_csv(f"{base}.CFPTMatrix.tsv", sep="\t", index_col=0)
    assert not np.allclose(cfpt.to_numpy(), cfpt.to_numpy().T)
    net = ResistanceNetworkSAMC(
        network=os.path.join(data, "nmt.network"), shapefile=None,
        sizes=None, coords=f"{base}.coords", variables=["LENGTH_KM"],
        inmat=f"{base}.ResistanceMatrix.tsv",
        agg_opts={"LENGTH_KM": "ARITH"}, out=str(tmp_path / "fit"),
        verbose=False)
    rate, res = net._fit_flow_alone()
    assert rate == pytest.approx(0.05, rel=1e-3)
    assert res["r2m"].iloc[0] > 0.999


def test_simulate_samc_matches_fitted_model(tmp_path):
    # The simulated distances are those of the same model as fitted: the
    # mean passage times of both directions, ordered by sample
    from resistnet.CFPT import CFPT
    spec = pd.DataFrame({"VAR": ["ele_mt_cav", "sgr_dk_rav"],
                         "WEIGHT": [1.0, 0.5], "TRANSFORM": [0, 5],
                         "SHAPE": [1, 3]})
    data, base = _simulate_samc(tmp_path, spec, flow=-0.02)
    sim = pd.read_csv(f"{base}.ResistanceMatrix.tsv", sep="\t", index_col=0)
    net = ResistanceNetworkSAMC(
        network=os.path.join(data, "nmt.network"), shapefile=None,
        sizes=None, coords=f"{base}.coords", variables=list(spec["VAR"]),
        inmat=None, agg_opts={v: "ARITH" for v in spec["VAR"]},
        out=str(tmp_path / "fit"), verbose=False)
    ind = []
    for v in net._predictors.columns:
        s = spec.set_index("VAR").loc[v]
        ind += [1, s["WEIGHT"], s["TRANSFORM"], s["SHAPE"], 0]
    T = CFPT(net._transition_matrix(ind + [20.0, -0.02]),
             net._site_states, net._R_lat)
    fitted = pd.DataFrame((T + T.T) / 2,
                          index=list(net._points_labels.values()),
                          columns=list(net._points_labels.values()))
    fitted.index = fitted.index.astype(str)
    fitted.columns = fitted.columns.astype(str)
    sim.index = sim.index.astype(str)
    sim.columns = sim.columns.astype(str)
    fitted = fitted.loc[sim.index, sim.columns]
    np.testing.assert_allclose(sim.to_numpy(),
                               fitted.to_numpy() / fitted.to_numpy().max(),
                               rtol=1e-12)


def test_ensemble_samc(samc_net, tmp_path):
    # ensembleResistnet.py averages SAMC runs as a SAMC model, and refuses
    # to mix them with standard runs
    import json
    import shutil
    import subprocess
    import sys
    from resistnet.CFPT import CFPT
    reps = tmp_path / "reps"
    reps.mkdir()
    for seed in (1, 2):
        runner = ModelRunnerTPE(samc_net, seed=seed, verbose=False)
        runner.run_tpe(max_evals=3, threads=1, reps=2, max_shape=10,
                       n_startup=2, out=str(reps / f"samc_{seed}"),
                       verbose=False, plot=False)
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                          "ensembleResistnet.py")
    args = ["-n", os.path.join(data, "nmt.network"), "-c", "HYRIV_ID",
            "-C", os.path.join(data, "nmt.pointCoords.txt"), "-t", "1",
            "-X"]
    out = str(tmp_path / "ens")
    p = subprocess.run([sys.executable, script, "-i", str(reps), "-o", out,
                        *args], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    with open(f"{out}.modelSpec.json") as fh:
        assert all(m["type"] == "samc" for m in json.load(fh)["models"])
    # The best model's passage times, as computed directly
    row = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t",
                      na_values="-").iloc[0]
    ind = []
    for v in samc_net._predictors.columns:
        ind += ([1, row[f"{v}_weight"], row[f"{v}_trans"],
                 row[f"{v}_shape"], row[f"{v}_asym"]] if row[v] == 1
                else [0, 1.0, 0, 1.0, 0])
    ind += [row["max_resistance"], row["flow"]]
    T = CFPT(samc_net._transition_matrix(ind),
             samc_net._site_states, samc_net._R_lat)
    written = pd.read_csv(f"{out}.Model-Best.CFPTMatrix.tsv", sep="\t",
                          index_col=0)
    np.testing.assert_allclose(written.to_numpy(), T, rtol=1e-10)

    # Mixed with a standard run: refused
    mixed = tmp_path / "mixed"
    mixed.mkdir()
    shutil.copy(reps / "samc_1.HallOfFame.tsv", mixed)
    shutil.copy(os.path.join(data, "nmt_ensemble", "replicates",
                             "nmt_1.HallOfFame.tsv"), mixed)
    p = subprocess.run([sys.executable, script, "-i", str(mixed), "-o",
                        str(tmp_path / "mix"), *args],
                       capture_output=True, text=True)
    assert p.returncode != 0
    assert "mix standard and SAMC" in p.stderr + p.stdout


def test_flow_profile(samc_net, tmp_path):
    # Flow alone at every rate tried, both directions: its best is the
    # samc_flow reference, and rate 0 is samc_distance
    out = str(tmp_path / "ref")
    table = samc_net.evaluate_null_model(out)
    profile = pd.read_csv(f"{out}.flowProfile.tsv", sep="\t")
    assert (profile["rate"] < 0).any() and (profile["rate"] > 0).any()
    assert profile["aic"].min() == pytest.approx(
        table.loc["samc_flow", "aic"], abs=1e-6)
    assert profile.loc[profile["rate"] == 0, "aic"].iloc[0] == \
        pytest.approx(table.loc["samc_distance", "aic"], abs=1e-6)
    assert profile["delta_aic"].min() == 0
