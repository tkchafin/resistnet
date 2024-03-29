import numpy as np
import pandas as pd
import rpy2.robjects as ro
import resistnet.utils as utils
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


def MLPE_R(X, Y, scale=True):
    """
    Perform Mixed Linear Model Pedigree Estimate (MLPE) using R libraries.

    Args:
        X: The first matrix for MLPE analysis.
        Y: The second matrix for MLPE analysis.
        scale (bool): If True, scales the Y matrix.

    Returns:
        The result of the MLPE analysis.
    """
    x = utils.get_lower_tri(X)
    y = utils.get_lower_tri(Y)

    # Make ID table
    npop = X.shape[0]
    ID = utils.to_from_(npop)

    # Make ZZ matrix
    ZZ = pd.DataFrame(ZZ_mat_(npop, ID))

    # Center and scale y
    if scale:
        y = (y - np.mean(y)) / np.std(y)

    # Add x and y to data
    ID["x"] = x
    ID["y"] = y

    stuff = """
    library(Matrix)
    library(lme4)
    library(MuMIn)
    MLPE <- function(ID, ZZ, REML=FALSE) {
        ZZ <- Matrix(data.matrix(ZZ), sparse=TRUE)
        mod_list <- list('full'=lme4:::lFormula(y ~ x + (1 | pop1),
        data = ID, REML = FALSE),
                         "null"=lme4:::lFormula(y ~ 1 + (1 | pop1),
                         data = ID, REML = FALSE))
        MODs <- lapply(mod_list, function(x) {
            x$reTrms$Zt <- ZZ
            dfun <- do.call(lme4:::mkLmerDevfun, x)
            opt <- lme4:::optimizeLmer(dfun)
            MOD <- (lme4:::mkMerMod(environment(dfun),
            opt, x$reTrms, fr = x$fr))
        })
        full <- summary(MODs$full)
        null <- summary(MODs$null)
        loglik <- c(full$logLik[[1]])
        null_loglik <- c(null$logLik[[1]])
        r2 <- c(MuMIn:::r.squaredGLMM(MODs$full)[[1]])
        aic <- c(full$AICtab[[1]])
        null_aic <- c(null$AICtab[[1]])
        deltaAIC <- c(null$AICtab[[1]] - full$AICtab[[1]])
        df <- data.frame(loglik, r2, -1*aic, deltaAIC,
        null_loglik, -1*null_aic)
        colnames(df) <- c("loglik", "r2m", "aic",
        "delta_aic_null", "loglik_null", "aic_null")
        return(df)
    }
    """

    r = ro.r
    r['options'](warn=-1)
    r(stuff)
    mlpe = ro.globalenv['MLPE']
    with localconverter(ro.default_converter + pandas2ri.converter):
        mlpe_res = mlpe(ID, ZZ)

    return mlpe_res


def ZZ_mat_(pops, id):
    """
    Create a matrix for ZZ calculations.

    Args:
        pops: The number of populations.
        id: The ID dataframe.

    Returns:
        The calculated ZZ matrix.
    """
    zz = np.zeros((pops, id.shape[0]))
    for i in range(pops):
        for j in range(id.shape[0]):
            # If ID row j contains pop i+1, set zz[j, i] to 1
            if i + 1 in list(id.iloc[j]):
                zz[i, j] = 1
    return zz

# def testSM():
# 	x = np.array([2.6407333, 0.6583007, 1.9629121, 0.8529997, 2.2592001,
# 2.9629032, 2.0796441, 2.4179196, 0.2154603, 2.5016938])
# 	y = np.array([3.6407333, 1.6583007, 1.5629121, 0.4529997,
# 2.0592001, 2.0629032, 2.9796441, 3.1179196, 1.2154603, 1.5016938])

# 	#make ID table
# 	ID = to_from_(5) #nrow(matrix)
# 	print(ID)

# 	#make ZZ matrix
# 	ZZ = ZZ_mat_(5, ID)
# 	print(ZZ)

# 	#center and scale x and y
# 	x = (x-np.mean(x))/np.std(x)
# 	y = (y-np.mean(y))/np.std(y)

# 	#add x and y to data
# 	ID["x"] = x
# 	ID["y"] = y
# 	print(ID)

# 	#get matrix describing random effects structure
# 	vcs = getVCM(ID)

# 	#fit mixed model
# 	#equivalent to lme4 y ~ x + (1|pop1)
# 	#how do incorporate the 'ZZ' part??
# 	model = smf.mixedlm("y ~ x", data=ID, groups="pop1")
# 	#model = mlm.MixedLM(endog=ID["y"].to_numpy(), exog=ID["x"].to_numpy(),
# groups=ID["pop1"].to_numpy(), exog_re=ZZ)
# 	#model = mlm.MixedLM(endog=ID["y"].to_numpy(), exog=ID["x"].to_numpy(),
# groups=ID["pop1"].to_numpy(), exog_vc=vcs)
# 	model_fit = model.fit()
# 	print(model)
# 	print(model_fit)
# 	print(model_fit.summary())

# def getVCM(df):
# 	def f(x):
# 		n = x.shape[0]
# 		g2 = x["pop2"]
# 		u = g2.unique()
# 		u.sort()
# 		uv = {v: k for k, v in enumerate(u)}
# 		mat = np.zeros((n, len(u)))
# 		for i in range(n):
# 		    mat[i, uv[g2[i]]] = 1
# 	    colnames = ["%d" % z for z in u]
# 		print(mat)
# 		return(mat, colnames)
# 	vcm = df.groupby("pop1").apply(f).to_list()
# 	print(vcm)
# 	mats = [x[0] for x in vcm]
# 	print(mats)
# 	colnames = [x[1] for x in vcm]
# 	names = ["pop2"]
# 	vcs = VCSpec(names, [colnames], [mats])
# 	return(vcs)
