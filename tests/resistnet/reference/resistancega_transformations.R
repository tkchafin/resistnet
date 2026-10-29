# ResistanceGA's transformations (vector versions, from its
# R/internal_helper_functions.R), applied to a test covariate: rescaled to
# 0-10 with SCALE.vector, transformed with parm = c(eq, shape, max = 1),
# then rescaled to 0-1 for comparison.
args <- commandArgs(trailingOnly = TRUE)
src <- readLines(args[1])
# Only the function definitions needed (no package dependencies)
defs <- c("SCALE.vector", "Monomolecular", "Inv.Monomolecular",
          "Inv.Rev.Monomolecular", "Rev.Monomolecular", "Ricker",
          "Inv.Ricker", "Inv.Rev.Ricker", "Rev.Ricker")
env <- new.env()
exprs <- parse(text = src)
for (e in exprs) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% defs) eval(e, env)
}
raw <- c(0, 0.3, 1.1, 2.5, 4, 7.25, 12, 15.5, 22, 30, 41.8, 55, 70, 88.2, 100)
x <- env$SCALE.vector(raw, 0, 10)
rescale01 <- function(v) (v - min(v)) / (max(v) - min(v))
out <- list()
for (f in setdiff(defs, "SCALE.vector")) {
  for (shape in c(0.5, 1, 3, 10, 50)) {
    y <- env[[f]](x, parm = c(1, shape, 1))
    out[[length(out) + 1]] <- data.frame(transformation = f, shape = shape,
                                         index = seq_along(x) - 1,
                                         raw = raw, value = rescale01(y))
  }
}
res <- do.call(rbind, out)
write.table(res, args[2], sep = "\t", row.names = FALSE, quote = FALSE)
cat(nrow(res), "rows\n")
