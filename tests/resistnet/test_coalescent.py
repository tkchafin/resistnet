import os

import numpy as np
import pandas as pd
import pytest

import resistnet.coalescent as co

REFERENCE = os.path.join(os.path.dirname(__file__), "reference")

# A small river (a tree of six demes: 0-1-2-3 main stem, 4-5 a tributary
# joining at 1), with a jump rate in each direction of every edge
# (asymmetric) and a coalescence rate per deme
TREE = [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5)]
RATES = {(0, 1): 0.7, (1, 0): 1.3, (1, 2): 0.4, (2, 1): 2.1, (2, 3): 1.0,
         (3, 2): 0.6, (1, 4): 0.9, (4, 1): 0.3, (4, 5): 1.7, (5, 4): 0.5}
GAMMA = np.array([1.0, 0.8, 1.5, 0.6, 1.2, 0.9])


def _edges():
    return list(RATES), np.array(list(RATES.values()))


def _divergences():
    """Divergences to fit: the model's coalescence times, perturbed by a
    fixed pattern (so that the residuals are not zero)."""
    edges, rates = _edges()
    C = co.coalescence_times(co.generator(rates, edges, 6), GAMMA)
    k = np.add.outer(np.arange(6), np.arange(6))
    return C * (1 + 0.03 * np.sin(k + 1))


def test_coalescence_times_equation():
    # C solves G C + C G' - diag(gamma) o C = -1 1' (their equation 7)
    edges, rates = _edges()
    G = co.generator(rates, edges, 6).toarray()
    np.testing.assert_allclose(G.sum(axis=1), 0, atol=1e-15)
    C = co.coalescence_times(G, GAMMA)
    lhs = G @ C + C @ G.T - np.diag(GAMMA * np.diag(C))
    np.testing.assert_allclose(lhs, -1, rtol=1e-12)
    np.testing.assert_allclose(C, C.T)
    # One deme: 1 / gamma
    assert co.coalescence_times(np.zeros((1, 1)), [2.0])[0, 0] == \
        pytest.approx(0.5)


def test_upper_order():
    # As R's M[upper.tri(M, diag = TRUE)]: column by column
    M = np.arange(9).reshape(3, 3)
    assert list(co.upper(M)) == [0, 1, 4, 2, 5, 8]


@pytest.mark.skipif(
    not os.path.exists(os.path.join(REFERENCE, "lundgren_ralph_H.tsv")),
    reason="reference/lundgren_ralph.R not run yet (needs their code)")
def test_matches_lundgren_ralph_code():
    # Their own HCalc and H_and_lllh (gene.flow.inference, run in R; see
    # reference/lundgren_ralph.R) on the same graph, rates, coalescence
    # rates and divergences
    edges, rates = _edges()
    H_ref = pd.read_csv(os.path.join(REFERENCE, "lundgren_ralph_H.tsv"),
                        sep="\t", header=None).to_numpy()
    C = co.coalescence_times(co.generator(rates, edges, 6), GAMMA)
    np.testing.assert_allclose(C, H_ref, rtol=1e-12)
    D = pd.read_csv(os.path.join(REFERENCE, "lundgren_ralph_inputs.tsv"),
                    sep="\t", header=None).to_numpy()
    np.testing.assert_allclose(D, _divergences(), rtol=1e-12)
    lllh = float(open(os.path.join(REFERENCE,
                                   "lundgren_ralph_lllh.txt")).read())
    lp = co.log_posterior(co.upper(D), co.upper(C), (D.mean() / 50) ** 2,
                          rates, GAMMA)
    assert lp == pytest.approx(lllh, rel=1e-12)
    model = co.Model(edges, 6, const_coal=False)
    np.testing.assert_allclose(model.hcalc(rates, GAMMA), co.upper(C))


def test_divergence():
    # Between groups: the share of loci at which two haplotypes differ,
    # averaged over all pairs; within: over pairs of distinct haplotypes
    rng = np.random.default_rng(1)
    X = (rng.random((50, 7)) < 0.4).astype(int)
    groups = np.array([0, 0, 0, 1, 1, 2, 2])
    D = co.divergence(X, groups)
    for a in range(3):
        for b in range(3):
            ia = np.flatnonzero(groups == a)
            ib = np.flatnonzero(groups == b)
            pairs = [(i, j) for i in ia for j in ib if i != j]
            ref = np.mean([np.mean(X[:, i] != X[:, j]) for i, j in pairs])
            assert D[a, b] == pytest.approx(ref)


def test_mh_recovers_rates():
    # Noise-free divergences of a two-deme graph with asymmetric rates:
    # the posterior concentrates near the truth (the start is the best of
    # prior draws; adaptive proposals)
    edges = [(0, 1), (1, 0)]
    rates = np.array([0.5, 2.0])
    C = co.coalescence_times(co.generator(rates, edges, 2), [1.0, 1.0])
    model = co.Model(edges, 2, const_coal=True)
    rng = np.random.default_rng(3)
    out = co.mh(model, C, iterations=6000, sig2=(C.mean() / 500) ** 2,
                rng=rng)
    g = out["g"][3000:].mean(axis=0)
    gam = out["gam"][3000:].mean()
    pred = co.coalescence_times(co.generator(g, edges, 2), [gam, gam])
    np.testing.assert_allclose(pred, C, rtol=0.02)


if __name__ == "__main__":
    # Writes the inputs read by reference/lundgren_ralph.R
    edges, rates = _edges()
    pd.DataFrame({"x": [e[0] for e in edges], "y": [e[1] for e in edges],
                  "rate": rates}).to_csv(
        os.path.join(REFERENCE, "lundgren_ralph_edges.tsv"), sep="\t",
        index=False)
    np.savetxt(os.path.join(REFERENCE, "lundgren_ralph_gamma.txt"), GAMMA,
               fmt="%.17g")
    np.savetxt(os.path.join(REFERENCE, "lundgren_ralph_inputs.tsv"),
               _divergences(), delimiter="\t", fmt="%.17g")


def _coal_net(tmp_path, cell_km, density=1000.0, variables=("dor_pc_pva",)):
    import resistnet
    from resistnet.coalescent_network import ResistanceNetworkCoalescent
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    v = list(variables)
    return ResistanceNetworkCoalescent(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"), sizes=None,
        coords=os.path.join(data, "nmt.pointCoords.txt"), variables=v,
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={x: "ARITH" for x in v},
        out=str(tmp_path / f"coal{cell_km:g}"), verbose=False,
        cell_km=cell_km, density=density)


def test_coalescent_network_symmetric_is_river_distance(tmp_path):
    # With symmetric movement (no flow), resistance 1 per km and high
    # gene flow, linearised FST from the coalescence times is river
    # distance up to a shift and scale (coalescence and commute times are
    # equal when hitting times are symmetric; Lundgren & Ralph 2019), so
    # MLPE fits it as well as the standard model's distance alone
    net = _coal_net(tmp_path, 2.0, density=1e5)
    lin = net.linfst(net._distance_rates(0.0))
    from resistnet.resistance_network import ResistanceNetwork
    import resistnet
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    std = ResistanceNetwork(
        network=None, shapefile=os.path.join(data, "nmt.gpkg"),
        coords=os.path.join(data, "nmt.pointCoords.txt"),
        variables=["dor_pc_pva"],
        inmat=os.path.join(data, "nmt.popGenDistMat.txt"),
        agg_opts={"dor_pc_pva": "ARITH"}, out=str(tmp_path / "std"),
        verbose=False)
    null = std.evaluate_null_model()
    res = net._flow_alone(0.0)
    assert float(res["loglik"].iloc[0]) == pytest.approx(
        float(null.loc["distance_only", "loglik"]), abs=0.05)
    assert lin.min() >= 0
    assert np.allclose(lin, lin.T)


def test_coalescent_network_spacing_and_direction(tmp_path):
    # The lattice's spacing barely changes the fit (finite-volume
    # movement and coalescence rates); a flow term changes it, and the
    # model's specification records the coalescent settings
    fits = {}
    for cell in (2.0, 5.0):
        net = _coal_net(tmp_path, cell)
        fits[cell] = [float(net._flow_alone(r)["loglik"].iloc[0])
                      for r in (0.0, -0.02)]
    np.testing.assert_allclose(fits[2.0], fits[5.0], atol=1.5)
    assert fits[2.0][1] != pytest.approx(fits[2.0][0])
    spec = net.model_spec([1, 1.0, 7, 1.0, 0, 50.0, -0.02])
    assert spec["type"] == "coalescent" and spec["density"] == 1000.0
    assert "absorption" not in spec


def test_run_coalescent_model(tmp_path):
    # runResistnet.py --model coalescent: its outputs, and the coalescent
    # reference models (distance alone, flow alone)
    import subprocess
    import sys
    import resistnet
    data = os.path.join(os.path.dirname(resistnet.__file__), "data")
    script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                          "runResistnet.py")
    out = str(tmp_path / "run" / "nmt")
    p = subprocess.run(
        [sys.executable, script, "-n", os.path.join(data, "nmt.network"),
         "-c", os.path.join(data, "nmt.pointCoords.txt"),
         "-g", os.path.join(data, "nmt.popGenDistMat.txt"),
         "--vars", "dor_pc_pva", "--model", "coalescent", "-i", "3", "-r",
         "1", "-S", "2", "-t", "1", "--seed", "1", "-o", out],
        capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    prefix = os.path.join(out, "coalescent", "nmt")
    for ext in ("HallOfFame.tsv", "Model-Best.CoalescenceMatrix.tsv",
                "Model-Average.TransitionEdges.tsv", "modelSpec.json",
                "pairwise.tsv"):
        assert os.path.exists(f"{prefix}.{ext}"), ext
    null = pd.read_csv(f"{prefix}.Null-Model.tsv", sep="\t", index_col=0)
    assert {"coal_distance", "coal_flow"} <= set(null.index)
    pw = pd.read_csv(f"{prefix}.pairwise.tsv", sep="\t")
    assert "best_model_distance" in pw.columns
