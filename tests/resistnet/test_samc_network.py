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
# each site to one site, with movement weights from a directional variable
# under the original (unbounded) resistance-to-conductance mapping. Plain
# sparse LU solves were off by up to 230% (runoff, to site 9) and by a
# factor of 1e31 (temperature, to site 5)
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
    """spec: {variable: (weight, transformation, shape, directional)}"""
    ind = []
    for v in VARIABLES:
        ind += [1, *spec[v]] if v in spec else [0, 1.0, 0, 1, 0]
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
    idx = {e: i for i, e in enumerate(net._edge_order)}
    root = sorted(idx[d[net.reachid_col]] for _, _, d in
                  net._K.in_edges(net._origin, data=True))
    assert sorted(net._edge_order[i] for i in root) == [70648307, 70648308]
    assert net._adj[root[0], root[1]] == -1
    assert net._adj[root[1], root[0]] == -1


def test_composite_surface_order_independent(samc_net):
    # A directional variable after a symmetric one used to lose its
    # orientation (CSR/CSC data arrays were added together)
    spec = {"run_mm_cyr": (0.6, 0, 1, 0), "sgr_dk_rav": (0.8, 0, 1, 1)}
    W = samc_net._build_composite_surface(model(spec)).toarray()

    perm = list(reversed(VARIABLES))
    rev = copy.copy(samc_net)
    rev._predictors = samc_net._predictors[perm]
    ind = []
    for v in perm:
        ind += [1, *spec[v]] if v in spec else [0, 1.0, 0, 1, 0]
    W_rev = rev._build_composite_surface(ind).toarray()
    np.testing.assert_array_equal(W, W_rev)


def test_directional_orientation(samc_net):
    # Moving towards a reach with a higher value is harder
    W = samc_net._build_composite_surface(
        model({"sgr_dk_rav": (1.0, 0, 1, 1)})).tocsr()
    v = samc_net._predictors["sgr_dk_rav"].to_numpy()
    rows, cols = samc_net._adj.nonzero()
    up = v[cols] > v[rows]
    w_ij = np.asarray(W[rows, cols]).ravel()
    w_ji = np.asarray(W[cols, rows]).ravel()
    assert up.any()
    assert np.all(w_ij[up] <= w_ji[up])
    assert np.mean(w_ij[up] < w_ji[up]) > 0.9


@pytest.mark.parametrize("var,site,ref", [
    ("run_mm_cyr", 9, REF_RUN_DIR_TO_SITE9),
    ("tmp_dc_cyr", 5, REF_TMP_DIR_TO_SITE5)])
def test_cfpt_precision(samc_net, var, site, ref):
    from resistnet.CFPT import CFPT, transition_matrix
    from resistnet.utils import minmax, minmax_nonzero
    v = samc_net._predictors[var].to_numpy(dtype=float)
    v = (v - v.min()) / (v.max() - v.min())
    rows, cols = samc_net._adj.nonzero()
    w = minmax_nonzero(1 / minmax_nonzero(minmax(v[cols] - v[rows])))
    W = sp.csr_matrix((w, (rows, cols)), shape=samc_net._adj.shape)
    T = CFPT(transition_matrix(W, samc_net._R), samc_net._edge_site_indices)
    np.testing.assert_allclose(T[:, site], ref, rtol=1e-10)


def test_resistance_mapping(samc_net):
    # Resistances run from 1 to the maximum resistance, so conductances
    # span exactly that contrast, and two-level variables are kept
    for var in ["run_mm_cyr", "dor_pc_pva"]:
        for k in [1.0, 37.0, 500.0]:
            W = samc_net._build_composite_surface(
                model({var: (1.0, 0, 1, 0)}) + [k])
            assert W.data.max() == pytest.approx(1.0)
            assert W.data.max() / W.data.min() == pytest.approx(k)
    two = copy.copy(samc_net)
    two._predictors = samc_net._predictors.copy()
    dor = two._predictors["dor_pc_pva"]
    two._predictors["dor_pc_pva"] = np.where(dor > dor.median(), 10.0, 0.0)
    W = two._build_composite_surface(
        model({"dor_pc_pva": (1.0, 0, 1, 0)}) + [100.0])
    assert len(np.unique(W.data)) == 3   # mean of two levels: 0, 0.5, 1
    # Without a maximum resistance in the model, the default is used
    assert samc_net.max_resistance(model({"run_mm_cyr": (1, 0, 1, 0)})) \
        == 100.0


def test_evaluate_and_output(samc_net):
    ind = model({"run_mm_cyr": (0.7, 3, 10, 0), "sgr_dk_rav": (0.4, 0, 1, 1)})
    fitness, res = samc_net.evaluate(ind)
    assert np.isfinite(fitness)
    # loglik, r2m, aic, delta_aic_null (as for the base model)
    assert np.all(np.isfinite(res[:4]))
    assert fitness == res[2]

    cfpt, trans = samc_net.model_output(ind)
    n = len(samc_net._edge_site_indices)
    assert cfpt.shape == (n, n)
    assert np.all(np.diag(cfpt) == 0)
    assert np.all(cfpt[~np.eye(n, dtype=bool)] >= 1)
    # Probabilities of moving downstream + upstream + absorption sum to 1
    total = np.nansum(trans, axis=1) + samc_net._R.ravel()
    np.testing.assert_allclose(total, 1)


def test_search_space(samc_net):
    runner = ModelRunnerTPE(samc_net, seed=1, verbose=False)
    runner.set_tpe_parameters(max_shape=10)
    runner.initialize_space()
    from hyperopt.pyll.stochastic import sample
    rng = np.random.default_rng(1)
    draws = [sample(runner.space, rng=rng) for _ in range(300)]
    assert {d["run_mm_cyr_asym"] for d in draws} == {0, 1}
    assert {d["flow_dir"] for d in draws} == {-1, 0, 1}
    rates = [d["flow_rate"] for d in draws]
    assert min(rates) >= 1e-3 and max(rates) <= samc_net.max_flow
    runner.set_tpe_parameters(max_shape=10, fixAsym=True, fixFlow=True)
    runner.initialize_space()
    draws = [sample(runner.space, rng=rng) for _ in range(50)]
    assert {d["run_mm_cyr_asym"] for d in draws} == {0}
    assert {d["flow_dir"] for d in draws} == {0}


def test_surface_parameters(samc_net):
    ind = model({"run_mm_cyr": (0.7, 3, 10, 0), "sgr_dk_rav": (0.4, 0, 1, 1)})
    # 2 x (transformation, shape) + 1 relative weight + 2 directions
    assert samc_net.surface_parameters(ind) == 7
    # + the maximum resistance, when it is part of the model
    assert samc_net.surface_parameters(ind + [50.0]) == 8
    net = copy.copy(samc_net)
    net.fixAsym = True
    assert net.surface_parameters(ind) == 5
    net.fix_resist = 50.0
    assert net.surface_parameters(ind + [50.0]) == 5


def test_run_tpe_samc(samc_net, tmp_path):
    runner = ModelRunnerTPE(samc_net, seed=1, verbose=False)
    out = str(tmp_path / "run")
    runner.run_tpe(max_evals=3, threads=1, reps=2, max_shape=10,
                   n_startup=2, out=out, verbose=False, plot=False)
    hof = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t", na_values="-")
    assert all(f"{v}_asym" in hof.columns for v in VARIABLES)
    for v in VARIABLES:
        assert hof[f"{v}_asym"].dropna().isin([0, 1]).all()
    assert hof["max_resistance"].between(1, samc_net.max_resist).all()
    assert hof["flow"].abs().le(samc_net.max_flow).all()
    # AIC differences from the reference models (positive = better)
    null = pd.read_csv(f"{out}.Null-Model.tsv", sep="\t", index_col=0)
    for col, ref in [("delta_aic_distance", "samc_distance"),
                     ("delta_aic_flow", "samc_flow")]:
        np.testing.assert_allclose(
            hof[col], null.loc[ref, "aic"] - hof["aic"], rtol=1e-9)
    assert os.path.exists(f"{out}.Model-Average.CFPTMatrix.tsv")
    assert os.path.exists(f"{out}.Model-Average.TransitionEdges.tsv")
    # The models' specifications, best first, and the best model's
    # effective directions
    import json
    with open(f"{out}.modelSpec.json") as fh:
        spec = json.load(fh)
    assert [m["rank"] for m in spec["models"]] == list(
        range(1, len(hof) + 1))
    np.testing.assert_allclose([m["aic"] for m in spec["models"]],
                               hof["aic"])
    best = pd.read_csv(f"{out}.bestModel.tsv", sep="\t")
    assert "Effective_direction" in best.columns


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
        verbose=False, ne=50
    )
    R = net._R.ravel()
    np.testing.assert_allclose(R[net._edge_site_indices], 1 / (2 * 50))
    assert np.count_nonzero(R) == len(net._edge_site_indices)


def test_null_models(samc_net):
    table = samc_net.evaluate_null_model()
    assert list(table.index) == [
        "distance_only", "samc_distance", "samc_flow", "null"]
    # Equal movement weights = any model with a maximum resistance of 1
    ind = model({"run_mm_cyr": (1.0, 4, 20, 1)}) + [1.0]
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


def test_effective_directions(samc_net):
    # The direction setting applies to the transformed values: increasing
    # transformations keep it, decreasing ones reverse it, and a humped
    # Ricker (peak within the data: runoff spans the range) gives a mix
    v = "sgr_dk_rav"
    share = {}
    for t, s in [(0, 1), (5, 100), (6, 100), (7, 100), (8, 100)]:
        d = samc_net.effective_directions(model({v: (1.0, t, s, 1)}))
        share[(t, s)] = d[v]
    assert share[(0, 1)] == 1.0
    assert share[(5, 100)] == 1.0 and share[(8, 100)] == 1.0
    assert share[(6, 100)] == 0.0 and share[(7, 100)] == 0.0
    mixed = samc_net.effective_directions(
        model({"run_mm_cyr": (1.0, 1, 3, 1)}))["run_mm_cyr"]
    assert 0.05 < mixed < 0.95
    # Symmetric variables have no direction
    assert samc_net.effective_directions(model({v: (1.0, 0, 1, 0)})) == {}
    # Values read back from the Hall of Fame are numpy scalars
    ind = model({v: (1.0, np.int64(7), np.int64(100), np.int64(1))})
    assert samc_net.effective_directions(ind)[v] == 0.0


def _raw(net):
    """Raw values of the predictors, from their 0-10 values and ranges."""
    raw = net._predictors.copy()
    for v, (lo, hi) in net._predictor_range.items():
        raw[v] = lo + raw[v] / 10 * (hi - lo)
    return raw


@pytest.mark.parametrize("t,s", [(0, 1), (1, 3), (2, 5), (3, 4), (4, 6),
                                 (5, 2), (6, 8), (7, 3), (8, 2)])
def test_model_spec_reproduces_weights(samc_net, t, s):
    # Applying the saved specification to the training reaches' raw values
    # gives exactly the fitted model's movement weights
    import resistnet.predict as pr
    ind = model({"sgr_dk_rav": (0.6, t, s, 1),
                 "run_mm_cyr": (0.4, 3, 7, 0)}) + [37.0, -0.04]
    spec = samc_net.model_spec(ind)
    W = samc_net._build_composite_surface(ind).tocsr()
    rows, cols = samc_net._adj.nonzero()
    up = samc_net._flow_pos[cols] - samc_net._flow_pos[rows]
    w = pr.samc_weights(_raw(samc_net), rows, cols, up, spec)
    np.testing.assert_allclose(w, np.asarray(W[rows, cols]).ravel(),
                               rtol=1e-12)


def test_predict_reach_layer(samc_net, tmp_path):
    # From the HydroRIVERS layer itself (each state is one reach), for the
    # reaches of the fitted network, the predicted weights are those of the
    # fitted model for every move
    import geopandas as gpd
    import json
    import subprocess
    import sys
    import resistnet.predict as pr
    ind = model({"sgr_dk_rav": (0.6, 7, 3, 1),
                 "run_mm_cyr": (0.4, 1, 4, 0)}) + [20.0, 0.03]
    spec = {"reachid_col": "HYRIV_ID", "length_col": "LENGTH_KM",
            "models": [samc_net.model_spec(ind)]}
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    reaches = gpd.read_file(os.path.join(data, "nmt.gpkg"))
    reaches = reaches[reaches["HYRIV_ID"].isin(samc_net._edge_order)]
    layer = str(tmp_path / "reaches.gpkg")
    reaches.to_file(layer, driver="GPKG")
    table, moves, clamped = pr.predict_reaches(reaches, spec)
    assert all(v == 0 for v in clamped.values())
    W = samc_net._build_composite_surface(ind).tocsr()
    idx = {e: i for i, e in enumerate(samc_net._edge_order)}
    a = np.array([idx[e] for e in moves["from"]])
    b = np.array([idx[e] for e in moves["to"]])
    assert len(moves) == 2 * (len(samc_net._edge_order) - 1)
    np.testing.assert_allclose(moves["weight"], np.asarray(W[a, b]).ravel(),
                               rtol=1e-12)
    # The command-line script writes the same values
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
    assert os.path.exists(f"{out}.predicted.gpkg")
    assert os.path.exists(f"{out}.predictedMoves.tsv")


def test_predict_clamps_to_training_range(samc_net):
    import resistnet.predict as pr
    spec = samc_net.model_spec(
        model({"sgr_dk_rav": (1.0, 5, 3, 1)}) + [10.0, 0.0])
    var = spec["variables"]["sgr_dk_rav"]
    lo, hi = var["min"], var["max"]
    inside = pr.variable_values([lo, hi], var)
    outside = pr.variable_values([lo - 5, hi + 5], var)
    np.testing.assert_allclose(outside, inside)
    assert pr.out_of_range([lo - 1, lo, hi, hi + 1], var) == 0.5


def test_flow_positions(samc_net):
    # Differences between adjacent reaches are the distance moved upstream
    # between midpoints: half of each reach's length (LENGTH_KM)
    eid = samc_net.reachid_col
    att = {d[eid]: (d["DIST_DN_KM"], d["LENGTH_KM"], d["NEXT_DOWN"])
           for _, _, d in samc_net._K.edges(data=True)}
    dd, L, nd = np.array([att[e] for e in samc_net._edge_order]).T
    rows, cols = samc_net._adj.nonzero()
    pos = samc_net._flow_pos
    order = list(samc_net._edge_order)
    down = np.array([order[c] == nd[r] for r, c in zip(rows, cols)])
    up = np.array([order[r] == nd[c] for r, c in zip(rows, cols)])
    assert np.all(down | up)
    step = (L[rows] + L[cols]) / 2
    np.testing.assert_allclose(pos[cols] - pos[rows],
                               np.where(up, step, -step), atol=1e-9)
    # That is the change in DIST_DN_KM + LENGTH_KM / 2, up to HydroRIVERS'
    # own rounding: DIST_DN_KM differs from the next reach's DIST_DN_KM +
    # LENGTH_KM by up to 0.32 km in nmt.gpkg
    below = {e: i for i, e in enumerate(order)}
    gap = max(abs(dd[i] - dd[below[nd[i]]] - L[below[nd[i]]])
              for i in range(len(order)) if nd[i] in below)
    assert gap < 0.5
    mid = dd + L / 2
    np.testing.assert_allclose(pos[cols] - pos[rows], mid[cols] - mid[rows],
                               atol=gap + 1e-9)


def test_flow_term(samc_net):
    rows, cols = samc_net._adj.nonzero()
    up = samc_net._flow_pos[cols] - samc_net._flow_pos[rows]
    for rate in [0.2, -0.05]:
        # Flow alone: weights are exp(-rate * km moved upstream)
        W = samc_net._build_composite_surface(model({}) + [10.0, rate])
        np.testing.assert_allclose(np.asarray(W[rows, cols]).ravel(),
                                   np.exp(-rate * up), rtol=1e-12)
        # With variables, it multiplies their weights; it is per km, so
        # not affected by the maximum resistance
        spec = model({"run_mm_cyr": (1.0, 0, 1, 0)})
        W0 = samc_net._build_composite_surface(spec + [10.0, 0.0])
        Wf = samc_net._build_composite_surface(spec + [10.0, rate])
        np.testing.assert_allclose(
            np.asarray(Wf[rows, cols]).ravel(),
            np.asarray(W0[rows, cols]).ravel() * np.exp(-rate * up),
            rtol=1e-12)
    # Neither variables nor flow: no model
    assert samc_net._build_composite_surface(model({}) + [10.0, 0.0]) is None
    fitness, _ = samc_net.evaluate(model({}) + [10.0, 0.0])
    assert fitness == float("-inf")


def test_surface_parameters_flow(samc_net):
    ind = model({"run_mm_cyr": (0.7, 3, 10, 0), "sgr_dk_rav": (0.4, 0, 1, 1)})
    assert samc_net.surface_parameters(ind + [50.0, 0.0]) == 8
    assert samc_net.surface_parameters(ind + [50.0, -0.1]) == 9
    # Flow alone: only the rate (no maximum resistance without variables)
    assert samc_net.surface_parameters(model({}) + [50.0, 0.1]) == 1


REFERENCE = os.path.join(os.path.dirname(__file__), "reference")
# Models with references from samc::cond_passage (reference/README.md):
# {variable: (weight, transformation, shape, directional)}, maximum
# resistance and flow rate; None = equal movement weights
R_SAMC_MODELS = {
    "equal": None,
    "sgr_sym": ({"sgr_dk_rav": (1.0, 0, 1.0, 0)}, 50.0, 0.0),
    "ele_dir_flow": ({"ele_mt_cav": (1.0, 0, 1.0, 1)}, 20.0, 0.02),
}


@pytest.fixture(scope="module")
def ref_net(tmp_path_factory):
    # The variables of the R reference models
    out = tmp_path_factory.mktemp("ref")
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    v = ["sgr_dk_rav", "ele_mt_cav"]
    return ResistanceNetworkSAMC(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"),
        sizes=os.path.join(data, "nmt.popSizes.txt"),
        coords=os.path.join(data, "nmt.pointCoords.txt"), variables=v,
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={x: "ARITH" for x in v}, out=str(out / "ref"),
        verbose=False)


def _ref_transition_matrix(net, name):
    from resistnet.CFPT import transition_matrix
    if R_SAMC_MODELS[name] is None:
        return transition_matrix(abs(net._adj).astype(float), net._R)
    spec, K, rate = R_SAMC_MODELS[name]
    ind = []
    for v in net._predictors.columns:
        ind += [1, *spec[v]] if v in spec else [0, 1.0, 0, 1.0, 0]
    return net._transition_matrix(ind + [K, rate])


@pytest.mark.parametrize("name", sorted(R_SAMC_MODELS))
def test_cfpt_matches_r_samc_nmt(ref_net, name):
    # Conditional passage times among the NMT sites, as samc::cond_passage
    from resistnet.CFPT import CFPT
    ref = np.loadtxt(os.path.join(REFERENCE, f"nmt_samc_{name}.tsv"))
    Q = _ref_transition_matrix(ref_net, name)
    T = CFPT(Q, ref_net._edge_site_indices, ref_net._R)
    np.testing.assert_allclose(T, ref, rtol=1e-10, atol=1e-12)


def test_mlpe_reference_inputs_nmt(ref_net):
    # The genetic distances and the SAMC model distances behind the lme4
    # references (test_MLPE.py) are what ResistNet computes for NMT
    from resistnet.CFPT import CFPT
    Y = np.loadtxt(os.path.join(REFERENCE, "nmt_mlpe_Y.tsv"))
    # (input parsing can differ in the last digit between pandas versions)
    np.testing.assert_allclose(ref_net._gendist, Y, rtol=1e-12, atol=0)
    X = np.loadtxt(os.path.join(REFERENCE,
                                "nmt_mlpe_X_samc_ele_dir_flow.tsv"))
    T = CFPT(_ref_transition_matrix(ref_net, "ele_dir_flow"),
             ref_net._edge_site_indices, ref_net._R)
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
                         "SHAPE": [1, 3], "ASYM": [1, 0]})
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
        ind += [1, s["WEIGHT"], s["TRANSFORM"], s["SHAPE"], s["ASYM"]]
    T = CFPT(net._transition_matrix(ind + [20.0, -0.02]),
             net._edge_site_indices, net._R)
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
    best = pd.read_csv(f"{out}.bestModel.tsv", sep="\t")
    assert "Effective_direction" in best.columns
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
             samc_net._edge_site_indices, samc_net._R)
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
