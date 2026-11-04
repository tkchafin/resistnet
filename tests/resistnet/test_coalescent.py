# Lundgren & Ralph's coalescence times (resistnet.coalescent), and the
# coalescent model on the NMT network
import os

import networkx as nx
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


def _river():
    """A small river network directed downstream (nodes (row, col),
    lengths in km): a main stem to the outlet at (6, 15), with tributaries
    joining at (6, 5) and (6, 11), and a degree-2 node at (6, 8) splitting
    the stretch between them in two reaches."""
    reaches = {1: ((6, 0), (6, 5)), 2: ((0, 5), (6, 5)),
               31: ((6, 5), (6, 8)), 32: ((6, 8), (6, 11)),
               4: ((11, 11), (6, 11)), 5: ((6, 11), (6, 15))}
    K = nx.DiGraph()
    for rid, (u, v) in reaches.items():
        L = abs(v[0] - u[0]) + abs(v[1] - u[1])
        K.add_edge(u, v, HYRIV_ID=rid, LENGTH_KM=float(L))
    return K, list(reaches)


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


def test_matches_lundgren_ralph_code():
    # Their own HCalc (gene.flow.inference, run in R; see
    # reference/lundgren_ralph.R) on the same graph, rates and coalescence
    # rates
    edges, rates = _edges()
    H_ref = pd.read_csv(os.path.join(REFERENCE, "lundgren_ralph_H.tsv"),
                        sep="\t", header=None).to_numpy()
    C = co.coalescence_times(co.generator(rates, edges, 6), GAMMA)
    np.testing.assert_allclose(C, H_ref, rtol=1e-12)


def test_segment_lattice():
    # Points at confluences, tips, the outlet and kept nodes, and about
    # cell_km apart in between, whatever the reaches' boundaries; a link's
    # resistance is the reaches' resistances in series over the river it
    # covers
    K, order = _river()
    seg = co.SegmentLattice(K, order, "HYRIV_ID", "LENGTH_KM", cell_km=2.0,
                            keep_nodes=[(6, 0)])
    # (6, 8) is not a point; the stretch (6, 5) - (6, 11) has 3 links
    with pytest.raises(KeyError):
        seg.states([(6, 8)])
    assert seg.overlap.sum() / 2 == pytest.approx(5 + 6 + 6 + 5 + 4)
    np.testing.assert_allclose(np.asarray(seg.overlap.sum(axis=1)).ravel(),
                               seg.length)
    # Resistance 1 per km: conductance 1 / length; otherwise in series
    W = seg.weights().tocsr()
    np.testing.assert_allclose(np.asarray(W[seg.a, seg.b]).ravel(),
                               1 / seg.length)
    r = np.array([2.0, 7.0, 1.0, 4.0, 12.0, 3.0])
    Wr = seg.weights(r).tocsr()
    np.testing.assert_allclose(np.asarray(Wr[seg.a, seg.b]).ravel(),
                               1 / (seg.overlap @ r))
    # The middle link of that stretch covers 1 km of each new reach
    mid = [i for i in range(len(seg.a))
           if seg.overlap[i, 2] > 0 and seg.overlap[i, 3] > 0]
    assert mid
    np.testing.assert_allclose(seg.overlap[mid[0], 2], 1.0)
    np.testing.assert_allclose(seg.overlap[mid[0], 3], 1.0)
    np.testing.assert_allclose(seg.pos[seg.states([(0, 5)])], 16.0)
    # One outlet only
    K.add_edge((20, 20), (21, 21), HYRIV_ID=9, LENGTH_KM=1.0)
    with pytest.raises(ValueError):
        co.SegmentLattice(K, order + [9], "HYRIV_ID", "LENGTH_KM")


def test_river_coalescent_by_hand():
    # Two tributaries A and B (2 km) join at c, 4 km above the outlet o;
    # points 2 km apart: A, B, c, m (the middle of c - o) and o. River
    # around each point: A 1, B 1, c 3, m 2, o 1 km
    K = nx.DiGraph()
    K.add_edge("A", "c", HYRIV_ID=10, LENGTH_KM=2.0)
    K.add_edge("B", "c", HYRIV_ID=11, LENGTH_KM=2.0)
    K.add_edge("c", "o", HYRIV_ID=12, LENGTH_KM=4.0)
    rc = co.RiverCoalescent(K, [10, 11, 12], ["A", "B", "o"], "HYRIV_ID",
                            "LENGTH_KM", cell_km=2.0, density=50.0)
    lat = rc.lattice
    A, B, c, o = lat.states(["A", "B", "c", "o"])
    m = ({*range(lat.n)} - {A, B, c, o}).pop()
    assert rc.sites == [A, B, o]
    pos = {A: 6.0, B: 6.0, c: 4.0, m: 2.0, o: 0.0}
    np.testing.assert_allclose(lat.pos[list(pos)], list(pos.values()))
    around = {A: 1.0, B: 1.0, c: 3.0, m: 2.0, o: 1.0}
    np.testing.assert_allclose(rc.around[list(around)],
                               list(around.values()))
    r, theta = np.array([2.0, 4.0, 1.5]), -0.1
    # Jump rate x -> y: 1 / (resistance per km x km) x exp(-theta x km
    # upstream), / river around x
    links = {(A, c): 2.0 * 2, (B, c): 4.0 * 2, (c, m): 1.5 * 2,
             (m, o): 1.5 * 2}
    rates = {}
    for (x, y), R in links.items():
        for u, v in ((x, y), (y, x)):
            rates[(u, v)] = (np.exp(-theta * (pos[v] - pos[u])) / R
                             / around[u])
    got = dict(zip(zip(lat.a, lat.b), rc.jump_rates(r, theta)))
    assert got.keys() == rates.keys()
    for k in rates:
        assert got[k] == pytest.approx(rates[k], rel=1e-14)
    # Equation 7 solved directly (dense), coalescence rates
    # 1 / (density x river around)
    n = lat.n
    G = np.zeros((n, n))
    for (u, v), q in rates.items():
        G[u, v] = q
    G -= np.diag(G.sum(axis=1))
    gamma = np.array([1 / (50.0 * around[i]) for i in range(n)])
    eye = np.eye(n)
    A_ = np.kron(G, eye) + np.kron(eye, G) - np.diag(
        np.where(np.eye(n).ravel() == 1, np.repeat(gamma, n), 0.0))
    C = np.linalg.solve(A_, -np.ones(n * n)).reshape(n, n)
    np.testing.assert_allclose(rc.coalescence_times(r, theta), C,
                               rtol=1e-10)
    # Linearised FST among the sites, (C_ij - Cw) / Cw
    s = [A, B, o]
    lin = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                w = (C[s[i], s[i]] + C[s[j], s[j]]) / 2
                lin[i, j] = (C[s[i], s[j]] - w) / w
    np.testing.assert_allclose(rc.linfst(r, theta), lin, rtol=1e-10)


def test_river_coalescent_spacing_and_flow():
    # Finite volumes: the lattice's spacing barely changes the result;
    # resistance 1 per km is distance alone; a flow term changes it
    K, order = _river()
    sites = [(6, 0), (0, 5), (11, 11), (6, 15), (6, 5)]
    r = np.array([2.0, 7.0, 1.0, 4.0, 12.0, 3.0])
    lin = {}
    for cell in (0.5, 1.0):
        rc = co.RiverCoalescent(K, order, sites, "HYRIV_ID", "LENGTH_KM",
                                cell_km=cell, density=100.0)
        lin[cell] = rc.linfst(r)
    np.testing.assert_allclose(lin[0.5], lin[1.0], rtol=0.01)
    np.testing.assert_allclose(rc.linfst(), rc.linfst(np.ones(6)),
                               rtol=1e-12)
    assert np.allclose(lin[1.0], lin[1.0].T)
    assert np.all(np.diag(lin[1.0]) == 0) and lin[1.0].min() >= 0
    flow = rc.linfst(r, -0.05)
    assert np.allclose(flow, flow.T)
    assert not np.allclose(flow, lin[1.0], rtol=1e-3)
    with pytest.raises(ValueError):
        co.RiverCoalescent(K, order, sites, "HYRIV_ID", "LENGTH_KM",
                           density=0)


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
    lin = net._lattice.linfst()
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


if __name__ == "__main__":
    # Writes the inputs read by reference/lundgren_ralph.R
    edges, rates = _edges()
    pd.DataFrame({"x": [e[0] for e in edges], "y": [e[1] for e in edges],
                  "rate": rates}).to_csv(
        os.path.join(REFERENCE, "lundgren_ralph_edges.tsv"), sep="\t",
        index=False)
    np.savetxt(os.path.join(REFERENCE, "lundgren_ralph_gamma.txt"), GAMMA,
               fmt="%.17g")
