library(Matrix)
library(lme4)
library(MuMIn)

MLPE <- function(ID, ZZ, REML=FALSE) {
  #print(ZZ)
	ZZ<-Matrix(data.matrix(ZZ), sparse=TRUE)
	#print(ZZ)
	# prepare model
  mod_list <- list('full'=lme4:::lFormula(y ~ x + (1 | pop1), data = ID, REML = FALSE),
                   "null"=lme4:::lFormula(y ~ 1 + (1 | pop1), data = ID, REML = FALSE))
  MODs <- lapply(mod_list, function(x) {
    x$reTrms$Zt <- ZZ
    # fit model
    dfun <- do.call(lme4:::mkLmerDevfun, x)
    opt <- lme4:::optimizeLmer(dfun)
    MOD <- (lme4:::mkMerMod(environment(dfun), opt, x$reTrms, fr = x$fr))
    #print(MOD)
  })
	models<-c("full")
	full<-summary(MODs$full)
	null<-summary(MODs$null)
	loglik<-c(full$logLik[[1]])
	r2<-c(MuMIn:::r.squaredGLMM(MODs$full)[[1]])
	aic<-c(full$AICtab[[1]])
	deltaAIC<-c(null$AICtab[[1]]-full$AICtab[[1]])
	df<-data.frame(loglik, r2, aic, deltaAIC)
	colnames(df) <- c("loglik", "r2m", "aic", "deltaAIC")
  # return model
  return(df)
}

