import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest

import resistnet.CFPT as cfpt
from resistnet.lattice import Lattice

REFERENCE = os.path.join(os.path.dirname(__file__), "reference")

# A small river drawn one cell wide on a 12 x 16 raster (rows, columns from
# 0): a main stem along row 6 to the outlet at (6, 15), with tributaries
# joining at (6, 5) and (6, 11). Each reach runs from its upstream node to
# its downstream node, with a resistance per km; cells are 1 km.
NROW, NCOL = 12, 16
REACHES = {
    1: ((6, 0), (6, 5), 2.0),
    2: ((0, 5), (6, 5), 7.0),
    3: ((6, 5), (6, 11), 1.5),
    4: ((11, 11), (6, 11), 12.0),
    5: ((6, 11), (6, 15), 3.0),
}
# Sites (nodes) and their effective sizes (absorption 1 / (2 Ne))
SITES = {(6, 0): 50, (0, 5): 80, (11, 11): 30, (6, 15): 100, (6, 5): 60}
# A directional variable (each reach's value) and the strength of its
# potential (U = STRENGTH x its values at the states)
DIRECTIONAL = {1: 0.2, 2: 0.9, 3: 0.4, 4: 0.0, 5: 0.7}
STRENGTH = 1.5


def river():
    """The river as a network directed downstream (nodes (row, col)),
    its reach order, and each state's cell on the raster."""
    K = nx.DiGraph()
    for rid, (u, v, _) in REACHES.items():
        L = abs(v[0] - u[0]) + abs(v[1] - u[1])
        K.add_edge(u, v, HYRIV_ID=rid, LENGTH_KM=float(L))
    order = list(REACHES)
    lat = Lattice(K, order, "HYRIV_ID", "LENGTH_KM", cell_km=1.0)
    cells = {i: v for i, v in enumerate(lat.nodes)}
    for e, states in enumerate(lat.reach_states):
        u, v, _ = REACHES[order[e]]
        L = len(states) + 1
        for k, s in enumerate(states, start=1):
            cells[s] = (u[0] + k * (v[0] - u[0]) // L,
                        u[1] + k * (v[1] - u[1]) // L)
    return K, order, lat, cells


def raster_inputs():
    """samc's rasters for the river: resistance per km of each cell (each
    reach's own; a confluence cell, the mean of the reaches meeting there)
    and absorption at the sites; NA off the river."""
    K, order, lat, cells = river()
    r = lat.state_resistance([REACHES[e][2] for e in order])
    res = np.full((NROW, NCOL), np.nan)
    ab = np.full((NROW, NCOL), np.nan)
    for s, (i, j) in cells.items():
        res[i, j] = r[s]
        ab[i, j] = 1 / (2 * SITES[(i, j)]) if (i, j) in SITES else 0.0
    return res, ab


def potential_raster():
    """The directional potential of each cell (see DIRECTIONAL); NA off
    the river."""
    K, order, lat, cells = river()
    U = STRENGTH * lat.state_values([DIRECTIONAL[e] for e in order])
    out = np.full((NROW, NCOL), np.nan)
    for s, (i, j) in cells.items():
        out[i, j] = U[s]
    return out


def _write(x, name):
    pd.DataFrame(x).to_csv(os.path.join(REFERENCE, name), sep="\t",
                           header=False, index=False, na_rep="NA")


def test_lattice_structure():
    K, order, lat, cells = river()
    # A state per km of river plus one per node
    km = sum(abs(v[0] - u[0]) + abs(v[1] - u[1])
             for u, v, _ in REACHES.values())
    assert lat.n == len(lat.nodes) + km - len(REACHES)
    np.testing.assert_allclose(lat.length, 1.0)
    # Each state is on its own cell, and the cells are distinct
    assert len(set(cells.values())) == lat.n
    # Distance to the outlet
    np.testing.assert_allclose(lat.pos[lat.states([(6, 15)])], 0.0)
    np.testing.assert_allclose(lat.pos[lat.states([(0, 5)])], 6 + 6 + 4)
    # Confluences: the mean of the three reaches meeting there
    r = lat.state_resistance([REACHES[e][2] for e in order])
    np.testing.assert_allclose(r[lat.states([(6, 5)])], (2 + 7 + 1.5) / 3)
    np.testing.assert_allclose(r[lat.states([(6, 11)])],
                               (1.5 + 12 + 3) / 3)
    np.testing.assert_allclose(r[lat.states([(6, 0)])], 2.0)
    # Neighbours are 1 cell apart on the raster
    for x, y in zip(lat.a, lat.b):
        (i, j), (k, m) = cells[x], cells[y]
        assert abs(i - k) + abs(j - m) == 1


def test_raster_inputs_unchanged():
    # The rasters samc_raster.R and samc_raster_dir.R read are those of
    # this river
    res, ab = raster_inputs()
    for x, name in [(res, "samc_raster_resistance.tsv"),
                    (ab, "samc_raster_absorption.tsv"),
                    (potential_raster(), "samc_raster_potential.tsv")]:
        ref = pd.read_csv(os.path.join(REFERENCE, name), sep="\t",
                          header=None).to_numpy(dtype=float)
        np.testing.assert_array_equal(np.isnan(ref), np.isnan(x))
        np.testing.assert_allclose(ref[~np.isnan(ref)], x[~np.isnan(x)],
                                   rtol=1e-14, atol=0)


def test_lattice_matches_samc_raster():
    # samc's model built by samc itself from the rasters (resistance,
    # absorption, 4-neighbours, conductance 1/mean(x); samc 3.2.1, see
    # reference/samc_raster.R) equals the lattice's: the transient matrix,
    # and the conditional passage times among the sites
    K, order, lat, cells = river()
    W = lat.weights([REACHES[e][2] for e in order])
    R = np.zeros(lat.n)
    sites = lat.states(list(SITES))
    R[sites] = [1 / (2 * ne) for ne in SITES.values()]
    Q = cfpt.transition_matrix(W, R).toarray()

    ref = pd.read_csv(os.path.join(REFERENCE, "samc_raster_Q.tsv"),
                      sep="\t", header=None).to_numpy(dtype=float)
    where = {(int(i) - 1, int(j) - 1): k
             for k, (i, j) in enumerate(ref[:, :2])}
    perm = [where[cells[s]] for s in range(lat.n)]
    np.testing.assert_allclose(Q, ref[:, 2:][np.ix_(perm, perm)],
                               rtol=1e-12, atol=1e-15)

    T = cfpt.CFPT(Q, sites, R)
    T_ref = pd.read_csv(os.path.join(REFERENCE, "samc_raster_cfpt.tsv"),
                        sep="\t", header=None).to_numpy(dtype=float)
    off = ~np.eye(len(sites), dtype=bool)
    np.testing.assert_allclose(T[off], T_ref[off], rtol=1e-10)


def test_lattice_matches_samc_raster_directional():
    # samc's model with a transition function of the order of the two
    # cells, exp(-(x_to - x_from)), on a raster of the directional
    # potential (samc 3.2.1, see reference/samc_raster_dir.R), equals the
    # lattice's with that potential and resistance 1 per km
    K, order, lat, cells = river()
    U = STRENGTH * lat.state_values([DIRECTIONAL[e] for e in order])
    W = lat.weights(potential=U)
    R = np.zeros(lat.n)
    sites = lat.states(list(SITES))
    R[sites] = [1 / (2 * ne) for ne in SITES.values()]
    Q = cfpt.transition_matrix(W, R).toarray()

    ref = pd.read_csv(os.path.join(REFERENCE, "samc_raster_dir_Q.tsv"),
                      sep="\t", header=None).to_numpy(dtype=float)
    where = {(int(i) - 1, int(j) - 1): k
             for k, (i, j) in enumerate(ref[:, :2])}
    perm = [where[cells[s]] for s in range(lat.n)]
    np.testing.assert_allclose(Q, ref[:, 2:][np.ix_(perm, perm)],
                               rtol=1e-12, atol=1e-15)
    # Not symmetric: the potential makes the two directions differ
    assert not np.allclose(Q, Q.T)

    T = cfpt.CFPT(Q, sites, R)
    T_ref = pd.read_csv(
        os.path.join(REFERENCE, "samc_raster_dir_cfpt.tsv"), sep="\t",
        header=None).to_numpy(dtype=float)
    off = ~np.eye(len(sites), dtype=bool)
    np.testing.assert_allclose(T[off], T_ref[off], rtol=1e-10)


def test_state_values():
    # Each reach's value at its middle, each node's the mean over the
    # half-links around it (all 1 km here), linear in between
    K, order, lat, cells = river()
    X = [DIRECTIONAL[e] for e in order]
    x = lat.state_values(X)
    up, down = lat.states([(6, 5), (6, 11)])
    np.testing.assert_allclose(x[up], (0.2 + 0.9 + 0.4) / 3)
    np.testing.assert_allclose(x[down], (0.4 + 0.0 + 0.7) / 3)
    np.testing.assert_allclose(x[lat.states([(0, 5)])], 0.9)
    # Reach 3, (6, 5) to (6, 11): 5 inner states, 1 to 5 km from (6, 5)
    inner = lat.reach_states[order.index(3)]
    want = [x[up] + (0.4 - x[up]) * k / 3 for k in (1, 2, 3)]
    want += [0.4 + (x[down] - 0.4) * k / 3 for k in (1, 2)]
    np.testing.assert_allclose(x[inner], want)
    # Constant values stay constant
    np.testing.assert_allclose(lat.state_values(np.full(5, 0.3)), 0.3)


def test_potential_is_flow_term():
    # The flow term is the potential rate x distance to the outlet
    K, order, lat, cells = river()
    r = [REACHES[e][2] for e in order]
    for rate in (0.3, -0.2):
        np.testing.assert_allclose(
            lat.weights(r, potential=rate * lat.pos).toarray(),
            lat.weights(r, rate).toarray(), rtol=1e-12)
        np.testing.assert_allclose(lat.reach_moves(0.0, rate * lat.pos),
                                   lat.reach_moves(rate), rtol=1e-12)
    # A potential multiplies each move by exp(-(U_to - U_from))
    U = STRENGTH * lat.state_values([DIRECTIONAL[e] for e in order])
    W0 = lat.weights(r).toarray()
    Wu = lat.weights(r, potential=U).toarray()
    for x, y in zip(lat.a, lat.b):
        assert Wu[x, y] == pytest.approx(W0[x, y] * np.exp(-(U[y] - U[x])))


def test_reach_moves_potential():
    # With a potential (and flow), the moves inside a reach are the mean
    # over its inner states of those of the transition matrix (resistance
    # 1 per km: the neighbouring nodes do not differ)
    K, order, lat, cells = river()
    U = STRENGTH * lat.state_values([DIRECTIONAL[e] for e in order])
    rate = 0.1
    Q = cfpt.transition_matrix(lat.weights(rate=rate, potential=U),
                               np.zeros(lat.n)).tocsr()
    moves = lat.reach_moves(rate, U)
    for e, u, v, states, _ in lat.chains:
        chain = [u] + list(states) + [v]
        down = [Q[chain[k], chain[k + 1]]
                for k in range(1, len(chain) - 1)]
        assert moves[e, 0] == pytest.approx(np.mean(down), rel=1e-12)
        assert moves[e, 1] == pytest.approx(1 - np.mean(down), rel=1e-12)
    # The variable's changes make moves differ between reaches
    assert np.ptp(moves[:, 0]) > 0.05


def test_weights_flow_and_distance():
    # Distance alone: conductance 1 / link length; the flow term
    # multiplies each move by exp(-rate x km upstream)
    K, order, lat, cells = river()
    W = lat.weights().tocoo()
    np.testing.assert_allclose(W.data, 1.0)
    rate = 0.3
    Wf = lat.weights(rate=rate).toarray()
    W0 = lat.weights().toarray()
    for x, y, up in zip(lat.a, lat.b, lat.upstream):
        assert Wf[x, y] == pytest.approx(W0[x, y] * np.exp(-rate * up))
    # Moves inside a reach, away from its ends: up and down equally likely
    # without flow; with it, as computed from the transition matrix at
    # such a state
    np.testing.assert_allclose(lat.reach_moves(0.0), 0.5)
    Q = cfpt.transition_matrix(lat.weights([REACHES[e][2] for e in order],
                                           rate), np.zeros(lat.n)).tocsr()
    moves = lat.reach_moves(rate)
    for e, states in enumerate(lat.reach_states):
        s = states[1]                   # neighbours: states[0], states[2]
        assert Q[s, states[2]] == pytest.approx(moves[e, 0])
        assert Q[s, states[0]] == pytest.approx(moves[e, 1])


def test_multi_part_reach():
    # A reach split into two edges with the same ID (as autoStreamTree
    # does with multi-part reaches): each part is its own chain, with the
    # reach's value at its middle; the symmetric model is unchanged
    K, order, lat, cells = river()
    K2 = nx.DiGraph()
    for rid, (u, v, _) in REACHES.items():
        L = abs(v[0] - u[0]) + abs(v[1] - u[1])
        if rid == 3:
            K2.add_edge(u, (6, 8), HYRIV_ID=3, LENGTH_KM=3.0)
            K2.add_edge((6, 8), v, HYRIV_ID=3, LENGTH_KM=3.0)
        else:
            K2.add_edge(u, v, HYRIV_ID=rid, LENGTH_KM=float(L))
    lat2 = Lattice(K2, order, "HYRIV_ID", "LENGTH_KM", cell_km=1.0)
    assert lat2.n == lat.n and len(lat2.chains) == 6
    e3 = order.index(3)
    parts = [c for c in lat2.chains if c[0] == e3]
    assert [len(c[3]) for c in parts] == [2, 2]
    assert sorted(parts[0][3] + parts[1][3]) == sorted(lat2.reach_states[e3])
    # Same movement weights and passage times between the sites
    r = [REACHES[e][2] for e in order]
    sites = list(SITES)
    for rate in (0.0, 0.2):
        T = [cfpt.CFPT(cfpt.transition_matrix(x.weights(r, rate),
                                              np.zeros(x.n)),
                       x.states(sites)) for x in (lat, lat2)]
        np.testing.assert_allclose(T[0], T[1], rtol=1e-12)
    # Values: the split node takes the reach's value (both half-links are
    # reach 3), and each part's inner states lie between its ends and the
    # reach's value
    X = [DIRECTIONAL[e] for e in order]
    x = lat2.state_values(X)
    assert x[lat2.states([(6, 8)])[0]] == pytest.approx(0.4)
    for _, u, v, states, _ in parts:
        lo, hi = min(x[u], x[v], 0.4), max(x[u], x[v], 0.4)
        assert np.all((x[states] >= lo - 1e-12) & (x[states] <= hi + 1e-12))
    assert np.all(np.isfinite(lat2.reach_moves(0.1, lat2.state_values(X))))


def test_one_outlet_required():
    K = nx.DiGraph()
    K.add_edge((0, 0), (0, 1), HYRIV_ID=1, LENGTH_KM=1.0)
    K.add_edge((5, 0), (5, 1), HYRIV_ID=2, LENGTH_KM=1.0)
    with pytest.raises(ValueError, match="outlets"):
        Lattice(K, [1, 2], "HYRIV_ID", "LENGTH_KM")


if __name__ == "__main__":
    # Writes the rasters and sites read by reference/samc_raster.R and
    # samc_raster_dir.R
    res, ab = raster_inputs()
    _write(res, "samc_raster_resistance.tsv")
    _write(ab, "samc_raster_absorption.tsv")
    _write(potential_raster(), "samc_raster_potential.tsv")
    pd.DataFrame([(i + 1, j + 1) for i, j in SITES],
                 columns=["row", "col"]).to_csv(
        os.path.join(REFERENCE, "samc_raster_sites.tsv"), sep="\t",
        index=False)
