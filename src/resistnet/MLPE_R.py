# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter
# import rpy2.robjects as ro

# def MLPE_R(Y, X, scale=True):
#     """
#     Perform Mixed Linear Model Pedigree Estimate (MLPE) using R libraries.

#     Args:
#         X: The first matrix for MLPE analysis.
#         Y: The second matrix for MLPE analysis.
#         scale (bool): If True, scales the Y matrix.

#     Returns:
#         The result of the MLPE analysis.
#     """
#     x = utils.get_lower_tri(X)
#     y = utils.get_lower_tri(Y)

#     # Make ID table
#     npop = X.shape[0]
#     ID = utils.to_from_(npop)

#     # Make ZZ matrix
#     ZZ = pd.DataFrame(ZZ_mat_(npop, ID))

#     # Center and scale y
#     if scale:
#         x = (x - np.mean(x)) / np.std(x)

#     # Add x and y to data
#     ID["x"] = x
#     ID["y"] = y

#     stuff = """
#     library(Matrix)
#     library(lme4)
#     library(MuMIn)
#     MLPE <- function(ID, ZZ, REML=FALSE) {
#         ZZ <- Matrix(data.matrix(ZZ), sparse=TRUE)
#         mod_list <- list('full'=lme4:::lFormula(y ~ x + (1 | pop1),
#         data = ID, REML = FALSE),
#                          "null"=lme4:::lFormula(y ~ 1 + (1 | pop1),
#                          data = ID, REML = FALSE))
#         MODs <- lapply(mod_list, function(x) {
#             x$reTrms$Zt <- ZZ
#             dfun <- do.call(lme4:::mkLmerDevfun, x)
#             opt <- lme4:::optimizeLmer(dfun)
#             MOD <- (lme4:::mkMerMod(environment(dfun),
#             opt, x$reTrms, fr = x$fr))
#         })
#         full <- summary(MODs$full)
#         null <- summary(MODs$null)
#         loglik <- c(full$logLik[[1]])
#         null_loglik <- c(null$logLik[[1]])
#         r2 <- c(MuMIn:::r.squaredGLMM(MODs$full)[[1]])
#         aic <- c(full$AICtab[[1]])
#         null_aic <- c(null$AICtab[[1]])
#         deltaAIC <- c(null$AICtab[[1]] - full$AICtab[[1]])
#         df <- data.frame(loglik, r2, -1*aic, deltaAIC,
#         null_loglik, -1*null_aic)
#         colnames(df) <- c("loglik", "r2m", "aic",
#         "delta_aic_null", "loglik_null", "aic_null")
#         return(df)
#     }
#     """

#     r = ro.r
#     r['options'](warn=-1)
#     r(stuff)
#     mlpe = ro.globalenv['MLPE']
#     with localconverter(ro.default_converter + pandas2ri.converter):
#         mlpe_res = mlpe(ID, ZZ)

#     return mlpe_res

# def ZZ_mat_(pops, id):
#     """
#     Create a matrix for ZZ calculations.

#     Args:
#         pops: The number of populations.
#         id: The ID dataframe.

#     Returns:
#         The calculated ZZ matrix.
#     """
#     zz = np.zeros((pops, id.shape[0]))
#     for i in range(pops):
#         for j in range(id.shape[0]):
#             # If ID row j contains pop i+1, set zz[j, i] to 1
#             if i + 1 in list(id.iloc[j]):
#                 zz[i, j] = 1
#     return zz
