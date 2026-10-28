import numpy as np
import pandas as pd
import pytest

import resistnet.MLPE as mlpe
import resistnet.utils as utils


def _mlpe_data(npop, seed=0):
    rng = np.random.default_rng(seed)
    pts = rng.random((npop, 2))
    X = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    noise = rng.random((npop, npop))
    Y = X + 0.3 * (noise + noise.T)
    np.fill_diagonal(Y, 0.0)
    return X, Y


def _cor_mlpe(npop):
    ids = utils.to_from_(npop)
    data = pd.DataFrame({"x": 0.0, "y": 0.0, "pop1": ids.pop1,
                         "pop2": ids.pop2})
    return mlpe.CorMLPE(form="~pop1 + pop2", data=data)


def test_lower_tri_matches_pair_ids():
    # each element of get_lower_tri must belong to the pair in to_from_
    for npop in (3, 4, 5, 12):
        mat = np.array(
            [[1000 * i + j for j in range(npop)] for i in range(npop)])
        lt = utils.get_lower_tri(mat)
        ids = utils.to_from_(npop)
        got = [{v // 1000, v % 1000} for v in lt]
        expected = [{a - 1, b - 1} for a, b in zip(ids.pop1, ids.pop2)]
        assert got == expected


def test_compute_ZZ():
    zz = _cor_mlpe(4)._compute_ZZ()
    ids = utils.to_from_(4)
    assert zz.shape == (4, 6)
    assert np.all(zz.sum(axis=0) == 2)
    for j, (a, b) in enumerate(zip(ids.pop1, ids.pop2)):
        assert zz[a - 1, j] == 1 and zz[b - 1, j] == 1


def test_compute_corr():
    cm = _cor_mlpe(5)
    corr = cm._compute_corr(cm._compute_ZZ(), 0.3)
    ids = utils.to_from_(5)
    pairs = list(zip(ids.pop1, ids.pop2))
    for i, p in enumerate(pairs):
        for j, q in enumerate(pairs):
            shared = len(set(p) & set(q))
            expected = 1.0 if i == j else (0.3 if shared == 1 else 0.0)
            assert corr[i, j] == pytest.approx(expected)


@pytest.mark.parametrize("rho", [0.001, 0.1, 0.25, 0.45])
def test_gls_mlpe_matches_statsmodels(rho):
    sm = pytest.importorskip("statsmodels.api")
    npop = 9
    X, Y = _mlpe_data(npop)
    x = utils.get_lower_tri(X)
    y = utils.get_lower_tri(Y)
    cm = _cor_mlpe(npop)
    ZZ = cm._compute_ZZ()
    corr = cm._compute_corr(ZZ, rho)

    design = np.column_stack((np.ones_like(x), x))
    ref = sm.GLS(y, design, sigma=corr).fit()
    beta, ssr, llf = mlpe.gls_mlpe(y, design, ZZ, rho)
    np.testing.assert_allclose(beta, ref.params, rtol=1e-8)
    assert ssr == pytest.approx(ref.ssr, rel=1e-8)
    assert llf == pytest.approx(ref.llf, rel=1e-8)

    ones = np.ones((len(y), 1))
    ref_null = sm.GLS(y, ones, sigma=corr).fit()
    _, _, llf_null = mlpe.gls_mlpe(y, ones, ZZ, rho)
    assert llf_null == pytest.approx(ref_null.llf, rel=1e-8)


# Reference values from ResistanceGA's MLPE.lmm (lme4, REML = FALSE), fit
# in R to Y = _mlpe_data(9, seed=0)[1] and the X given by each key:
# (loglik, loglik_null, AIC, MuMIn::r.squaredGLMM marginal R2)
LME4_REFERENCE = {
    "seed0": (28.98330005879670, -6.54558068896819, -49.9666001175934,
              0.8675290440366746),
    "seed1": (-6.25835158373607, -6.54558068896819, 20.5167031674721,
              0.0191385549049348),
    "sqrt_seed0": (26.25031359054580, -6.54558068896819, -44.5006271810917,
                   0.8457126480890063),
}


@pytest.mark.parametrize("case", sorted(LME4_REFERENCE))
def test_mlpe_matches_resistancega_lme4(case):
    X0, Y = _mlpe_data(9, seed=0)
    X = {"seed0": X0,
         "seed1": _mlpe_data(9, seed=1)[0],
         "sqrt_seed0": np.sqrt(X0)}[case]
    loglik, loglik_null, aic, r2m = LME4_REFERENCE[case]
    res = mlpe.MLPE(Y, X, scale=True).iloc[0]
    # seed0 and sqrt_seed0 have the population variance on the boundary
    # (rho = 0), where the 1e-8 lower bound on rho costs ~1e-7
    assert res["loglik"] == pytest.approx(loglik, abs=1e-6)
    assert res["loglik_null"] == pytest.approx(loglik_null, abs=1e-6)
    assert -res["aic"] == pytest.approx(aic, abs=1e-5)
    assert res["r2m"] == pytest.approx(r2m, abs=1e-6)


def test_mlpe_outputs():
    X, Y = _mlpe_data(10)
    res = mlpe.MLPE(Y, X, scale=True)
    assert list(res.columns) == [
        "loglik", "r2m", "aic", "delta_aic_null", "loglik_null", "aic_null"
    ]
    row = res.iloc[0]
    # aic is returned negated (higher is better); as lme4, k counts the
    # fixed effects plus two variance parameters
    assert -row["aic"] == pytest.approx(8 - 2 * row["loglik"])
    assert -row["aic_null"] == pytest.approx(6 - 2 * row["loglik_null"])
    assert row["delta_aic_null"] == pytest.approx(
        row["aic"] - row["aic_null"])
    assert 0.0 <= row["r2m"] <= 1.0
    # distance should explain the simulated data better than the null
    assert row["loglik"] > row["loglik_null"]
