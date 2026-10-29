# ResistanceGA's MLPE.lmm (its lower, To.From.ID and ZZ.mat, lFormula with
# the Zt swap, ML fit), fitted to the NMT genetic distances
# (nmt_mlpe_Y.tsv) and each model distance (nmt_mlpe_X_*.tsv), with the
# marginal R2 from MuMIn. Run in this directory: writes nmt_mlpe_lme4.tsv.
suppressPackageStartupMessages({library(lme4); library(Matrix)})

lower <- function(matrix) matrix[lower.tri(matrix)]
To.From.ID <- function(sampled_pops) {
  tmp <- matrix(0, nrow = sampled_pops, ncol = sampled_pops)
  colmat <- col(tmp); rowmat <- row(tmp)
  colmat[sampled_pops, 1] <- sampled_pops
  rowmat[sampled_pops, 1] <- 1
  data.frame(pop1 = factor(lower(colmat)), pop2 = factor(lower(rowmat)))
}
ZZ.mat <- function(ID) {
  Zl <- lapply(c("pop1", "pop2"), function(nm)
    Matrix::fac2sparse(ID[[nm]], "d", drop = FALSE))
  Reduce("+", Zl[-1], Zl[[1]])
}
fit <- function(form, dat, ZZ) {
  mod <- lme4::lFormula(form, data = dat, REML = FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(lme4::mkLmerDevfun, mod)
  opt <- lme4::optimizeLmer(dfun)
  lme4::mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr)
}

Y <- as.matrix(read.table("nmt_mlpe_Y.tsv", sep = "\t"))
ID <- To.From.ID(nrow(Y)); ZZ <- ZZ.mat(ID)
out <- list()
for (s in sub("^nmt_mlpe_X_(.*)\\.tsv$", "\\1",
              list.files(".", pattern = "^nmt_mlpe_X_.*\\.tsv$"))) {
  X <- as.matrix(read.table(paste0("nmt_mlpe_X_", s, ".tsv"), sep = "\t"))
  dat <- data.frame(ID, resistance = scale(lower(X)), response = lower(Y))
  full <- fit(response ~ resistance + (1 | pop1), dat, ZZ)
  null <- fit(response ~ 1 + (1 | pop1), dat, ZZ)
  out[[s]] <- data.frame(
    surface = s, loglik = as.numeric(logLik(full)),
    loglik_null = as.numeric(logLik(null)), AIC = AIC(full),
    r2m_MuMIn = MuMIn::r.squaredGLMM(full)[[1]])
}
write.table(do.call(rbind, out), "nmt_mlpe_lme4.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
