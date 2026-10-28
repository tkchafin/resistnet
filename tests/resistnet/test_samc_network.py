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
# the conditional passage times from each site to site 8, under the model
# with run_mm_cyr as a directional variable; without iterative refinement,
# the LU solves were off by up to 0.5%
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
    Q = samc_net._transition_matrix(model({"run_mm_cyr": (1.0, 0, 1, 1)}))
    from resistnet.CFPT import CFPT
    T = CFPT(Q, samc_net._edge_site_indices)
    np.testing.assert_allclose(T[:, 8], REF_RUN_DIR_TO_SITE8, rtol=1e-10)


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
    net = copy.copy(samc_net)
    net.fixAsym = True
    assert net.surface_parameters(ind) == 5


def test_run_tpe_samc(samc_net, tmp_path):
    runner = ModelRunnerTPE(samc_net, seed=1, verbose=False)
    out = str(tmp_path / "run")
    runner.run_tpe(max_evals=3, threads=1, reps=2, max_shape=10,
                   n_startup=2, out=out, verbose=False, plot=False)
    hof = pd.read_csv(f"{out}.HallOfFame.tsv", sep="\t")
    assert all(f"{v}_asym" in hof.columns for v in VARIABLES)
    assert os.path.exists(f"{out}.Model-Average.CFPTMatrix.tsv")
    assert os.path.exists(f"{out}.Model-Average.TransitionEdges.tsv")
