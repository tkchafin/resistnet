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

# 50-digit reference (sparse Gaussian elimination in decimal arithmetic) for
# the conditional passage times from each site to site 8, with movement
# weights from directional run_mm_cyr under the original (unbounded)
# resistance-to-conductance mapping; without iterative refinement, the LU
# solves were off by up to 0.5%
REF_RUN_DIR_TO_SITE8 = [
    4534.0497886998983, 4569.7948874062749, 4326.7003348945245,
    4118.7409360978836, 223.04489990490589, 180.64915345825153,
    117.7144752311454, 116.7144752311454, 0, 5350.7190134376415,
    6623.433419843991, 8000.5798515747247, 8224.3039948091573,
    8225.3039948091573, 9841.268509534475, 9498.6176551993849,
    9326.1171236755144, 8978.5337755095006, 8622.6503261687012,
    9264.5376148914402, 9259.6817896965149, 9230.3474544090604,
    9812.5401301094553, 8954.6404006587454, 8802.7583717583202,
    10007.857571594417, 9110.3217981500547, 8670.5813075881852,
    8377.8908993257264, 8340.6091647085377
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
        network=None, shapefile=os.path.join(data, "test.shp"),
        sizes=os.path.join(data, "test.popSizes.txt"),
        coords=os.path.join(data, "test.pointCoords.txt"),
        variables=VARIABLES,
        inmat=os.path.join(data, "test.popGenDistMat.txt"),
        agg_opts={v: "ARITH" for v in VARIABLES},
        out=str(out / "samc"), verbose=False
    )


def test_outlet_reaches_linked(samc_net):
    # test.shp is clipped at the outlet, so the reaches meeting there are
    # linked directly; otherwise the network splits in two
    n_comp, _ = sp.csgraph.connected_components(samc_net._adj,
                                                directed=False)
    assert n_comp == 1
    idx = {e: i for i, e in enumerate(samc_net._edge_order)}
    root = [idx[d[samc_net.reachid_col]] for _, _, d in
            samc_net._K.in_edges(samc_net._origin, data=True)]
    assert len(root) == 2
    assert samc_net._adj[root[0], root[1]] == -1
    assert samc_net._adj[root[1], root[0]] == -1


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


def test_cfpt_precision(samc_net):
    from resistnet.CFPT import CFPT, transition_matrix
    from resistnet.utils import minmax, minmax_nonzero
    v = samc_net._predictors["run_mm_cyr"].to_numpy(dtype=float)
    v = (v - v.min()) / (v.max() - v.min())
    rows, cols = samc_net._adj.nonzero()
    w = minmax_nonzero(1 / minmax_nonzero(minmax(v[cols] - v[rows])))
    W = sp.csr_matrix((w, (rows, cols)), shape=samc_net._adj.shape)
    T = CFPT(transition_matrix(W, samc_net._R), samc_net._edge_site_indices)
    np.testing.assert_allclose(T[:, 8], REF_RUN_DIR_TO_SITE8, rtol=1e-10)


def test_resistance_mapping(samc_net):
    # Resistances run from 1 to the maximum resistance, so conductances
    # span exactly that contrast, and two-level variables are kept
    for var in ["run_mm_cyr", "dor_pc_pva"]:
        for k in [1.0, 37.0, 500.0]:
            W = samc_net._build_composite_surface(
                model({var: (1.0, 0, 1, 0)}) + [k])
            assert W.data.max() == pytest.approx(1.0)
            assert W.data.max() / W.data.min() == pytest.approx(k)
    W = samc_net._build_composite_surface(
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
    hof = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t")
    assert all(f"{v}_asym" in hof.columns for v in VARIABLES)
    assert hof["max_resistance"].between(1, samc_net.max_resist).all()
    assert os.path.exists(f"{out}.Model-Average.CFPTMatrix.tsv")
    assert os.path.exists(f"{out}.Model-Average.TransitionEdges.tsv")


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
        network=None, shapefile=os.path.join(data, "test.shp"), sizes=None,
        coords=os.path.join(data, "test.pointCoords.txt"),
        variables=VARIABLES[:1], inmat=None,
        agg_opts={VARIABLES[0]: "ARITH"}, out=str(tmp_path / "ne"),
        verbose=False, ne=50
    )
    R = net._R.ravel()
    np.testing.assert_allclose(R[net._edge_site_indices], 1 / (2 * 50))
    assert np.count_nonzero(R) == len(net._edge_site_indices)
