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
    # The rasters samc_raster.R read are those of this river
    res, ab = raster_inputs()
    for x, name in [(res, "samc_raster_resistance.tsv"),
                    (ab, "samc_raster_absorption.tsv")]:
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


def test_one_outlet_required():
    K = nx.DiGraph()
    K.add_edge((0, 0), (0, 1), HYRIV_ID=1, LENGTH_KM=1.0)
    K.add_edge((5, 0), (5, 1), HYRIV_ID=2, LENGTH_KM=1.0)
    with pytest.raises(ValueError, match="outlets"):
        Lattice(K, [1, 2], "HYRIV_ID", "LENGTH_KM")


if __name__ == "__main__":
    # Writes the rasters and sites read by reference/samc_raster.R
    res, ab = raster_inputs()
    _write(res, "samc_raster_resistance.tsv")
    _write(ab, "samc_raster_absorption.tsv")
    pd.DataFrame([(i + 1, j + 1) for i, j in SITES],
                 columns=["row", "col"]).to_csv(
        os.path.join(REFERENCE, "samc_raster_sites.tsv"), sep="\t",
        index=False)
