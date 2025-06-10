import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize

import resistnet.utils as utils


# NOTE: This uses a fixed form now but will expand later
# form is not really used and is a placeholder
class CorMLPE:
    def __init__(self, form, data, rho=0.1):
        self.form = form
        self.data = data
        self.rho = rho

    def fit(self):
        x = self.data['x']
        y = self.data['y']
        ZZ = self._compute_ZZ()

        def objective(rho):
            corr = self._compute_corr(ZZ, rho)
            ols_model = sm.GLS(y, sm.add_constant(x), sigma=corr)
            result = ols_model.fit()
            return -result.llf

        bounds = [(0.0001, 0.4999)]
        result = minimize(objective, self.rho, bounds=bounds,
                        method='L-BFGS-B')
        self.rho = result.x[0]
        return self.rho

    def _compute_ZZ(self):
        pops = self.data['pop1'].nunique()
        id = self.data[['pop1', 'pop2']]
        zz = np.zeros((pops, id.shape[0]))
        for i in range(pops):
            for j in range(id.shape[0]):
                if i + 1 in list(id.iloc[j]):
                    zz[i, j] = 1
        return zz

    def _compute_corr(self, ZZ, rho):
        Z_scaled = ZZ * rho
        corr = Z_scaled.T @ Z_scaled
        np.fill_diagonal(corr, 1)
        return corr


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

    # Compute the correlation matrix
    ZZ = cor_mlpe._compute_ZZ()
    corr = cor_mlpe._compute_corr(ZZ, rho)

    # Fit the GLS model with the optimized correlation structure
    gls_model = sm.GLS(y, sm.add_constant(x), sigma=corr)
    result = gls_model.fit()

    # Fit the null model (intercept-only)
    null_model = sm.GLS(y, np.ones_like(y), sigma=corr)
    null_result = null_model.fit()

    # Calculate marginal R-squared
    var_fixed = np.var(result.fittedvalues)
    var_residual = result.mse_resid
    marginal_r2 = var_fixed / (var_fixed + var_residual)

    # Extract the desired statistics
    loglik = result.llf
    r2m = marginal_r2
    aic = result.aic
    loglik_null = null_result.llf
    aic_null = null_result.aic
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