import numpy as np
import pandas as pd
from scipy.optimize import minimize

import resistnet.utils as utils


# NOTE: This uses a fixed form now but will expand later
# form is not really used and is a placeholder
class CorMLPE:
    def __init__(self, form, data, rho=0.1):
        self.form = form
        self.data = data
        self.rho = rho

    def fit(self, X=None):
        # X: optional design matrix (default: intercept + x)
        y = self.data['y'].to_numpy(dtype=float)
        if X is None:
            x = self.data['x'].to_numpy(dtype=float)
            X = np.column_stack((np.ones_like(x), x))
        ZZ = self._compute_ZZ()
        self.ZZ = ZZ

        def objective(rho):
            return -gls_mlpe(y, X, ZZ, rho[0])[2]

        bounds = [(1e-8, 0.4999)]
        result = minimize(objective, self.rho, bounds=bounds,
                          method='L-BFGS-B')
        self.rho = result.x[0]
        return self.rho

    def _compute_ZZ(self):
        # pops x pairs incidence matrix: 1 where a pop is in a pair
        p1 = self.data['pop1'].to_numpy(dtype=int) - 1
        p2 = self.data['pop2'].to_numpy(dtype=int) - 1
        pops = max(p1.max(), p2.max()) + 1
        pairs = np.arange(len(p1))
        zz = np.zeros((pops, len(p1)))
        zz[p1, pairs] = 1
        zz[p2, pairs] = 1
        return zz

    def _compute_corr(self, ZZ, rho):
        # Dense pairs x pairs correlation: rho for pairs sharing a pop
        corr = rho * (ZZ.T @ ZZ)
        np.fill_diagonal(corr, 1)
        return corr


def gls_mlpe(y, X, ZZ, rho):
    """
    GLS fit under the MLPE correlation structure, without building the
    pairs x pairs correlation matrix.

    The correlation is Sigma = (1 - 2*rho) * I + rho * ZZ'ZZ (1 on the
    diagonal, rho for pairs sharing a population), so Sigma^-1 and
    log|Sigma| can be computed from the pops x pops matrix
    K = (1 - 2*rho) / rho * I + ZZ ZZ' (Woodbury identity and matrix
    determinant lemma). Results match statsmodels GLS with
    sigma=Sigma.

    Args:
        y (np.ndarray): Response (pairs,).
        X (np.ndarray): Design matrix (pairs, k).
        ZZ (np.ndarray): Pops x pairs incidence matrix.
        rho (float): Correlation for pairs sharing a population, in
                     (0, 0.5).

    Returns:
        tuple: (beta, ssr, llf), where ssr is the whitened residual sum of
               squares and llf the log-likelihood.
    """
    n = len(y)
    a = 1.0 - 2.0 * rho
    ZZt = ZZ @ ZZ.T
    K = (a / rho) * np.eye(ZZ.shape[0]) + ZZt

    def sigma_inv_inner(U, V):
        # U' Sigma^-1 V
        ZU = ZZ @ U
        ZV = ZZ @ V
        return (U.T @ V - ZU.T @ np.linalg.solve(K, ZV)) / a

    beta = np.linalg.solve(sigma_inv_inner(X, X), sigma_inv_inner(X, y))
    resid = y - X @ beta
    ssr = sigma_inv_inner(resid, resid)

    _, logdet = np.linalg.slogdet(np.eye(ZZ.shape[0]) + (rho / a) * ZZt)
    logdet += n * np.log(a)

    llf = -0.5 * n * (np.log(ssr) + 1 + np.log(2 * np.pi / n)) - 0.5 * logdet
    return beta, ssr, llf


def optimize_rho(data):
    cor_mlpe = CorMLPE(form='~pop1 + pop2', data=data)
    rho = cor_mlpe.fit()
    return rho, cor_mlpe


def MLPE(Y, X, scale=True):
    x = utils.get_lower_tri(X)
    y = utils.get_lower_tri(Y)

    # Make ID table
    npop = X.shape[0]
    ID = utils.to_from_(npop)

    # Center and scale x
    if scale:
        x = (x - np.mean(x)) / np.std(x)

    data = pd.DataFrame(
        {'x': x, 'y': y, 'pop1': ID['pop1'], 'pop2': ID['pop2']}
    )

    # Optimize rho and get the correlation structure
    rho, cor_mlpe = optimize_rho(data)

    # Fit the GLS model with the optimized correlation structure
    X_full = np.column_stack((np.ones_like(x), x))
    beta, ssr, loglik = gls_mlpe(y, X_full, cor_mlpe.ZZ, rho)

    # Fit the null model (intercept-only), with its own rho
    X_null = np.ones((len(y), 1))
    rho_null = CorMLPE(form='~pop1 + pop2', data=data).fit(X_null)
    _, _, loglik_null = gls_mlpe(y, X_null, cor_mlpe.ZZ, rho_null)

    # Calculate marginal R-squared as MuMIn::r.squaredGLMM (used by
    # ResistanceGA): var_f / (var_f + var_pop + var_resid), with the ML
    # variance components var_pop = rho * s2, var_resid = (1 - 2*rho) * s2
    var_fixed = np.var(X_full @ beta, ddof=1)
    s2 = ssr / len(y)
    var_pop = rho * s2
    var_residual = (1 - 2 * rho) * s2
    marginal_r2 = var_fixed / (var_fixed + var_pop + var_residual)

    # Extract the desired statistics
    # AIC counts the two variance parameters (rho and sigma), as in lme4
    r2m = marginal_r2
    aic = 2 * (X_full.shape[1] + 2) - 2 * loglik
    aic_null = 2 * (X_null.shape[1] + 2) - 2 * loglik_null
    delta_aic_null = aic_null - aic

    return pd.DataFrame({
        "loglik": [loglik],
        "r2m": [r2m],
        "aic": [-aic],
        "delta_aic_null": [delta_aic_null],
        "loglik_null": [loglik_null],
        "aic_null": [-aic_null]
    })

#     # Create the 'groups' column for the pairwise comparisons
#     ID['group'] = ID.apply(lambda row: f"{row['pop1']}", axis=1)

#     # Convert to DataFrame
#     data = pd.DataFrame({'x': x, 'y': y, 'group': ID['group']})

#     # # Fit mixed-effects model using Maximum Likelihood
#     # model = sm.MixedLM.from_formula("y ~ x", data, groups=data["group"],
#                exog_re=ZZ)
#     # result = model.fit(reml=False)

#     # # Fit the null model (intercept-only) using Maximum Likelihood
#     # null_model = sm.MixedLM.from_formula("y ~ 1", data,
#                       groups=data["group"], exog_re=ZZ)
#     # null_result = null_model.fit(reml=False)

#     # # Calculate marginal R-squared
#     # var_fixed = np.var(result.fittedvalues)
#     # var_random = result.cov_re.iloc[0, 0]
#     # var_residual = result.scale
#     # marginal_r2 = var_fixed / (var_fixed + var_random + var_residual)
