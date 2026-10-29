# ResistanceGA's composite resistance surface (as in its
# R/Combine_Surfaces.R) for the NMT reaches: each continuous surface
# rescaled to 0-10 with SCALE.vector and transformed with parm =
# c(eq, shape, max) (vector versions of its transformation functions, from
# R/internal_helper_functions.R), then summed and divided by its minimum
# ("rescale to min of 1"). Also each surface's range at max = 1, A (its
# values are 1 + max * A * f, f from 0 to 1).
# Usage: Rscript resistancega_surfaces.R internal_helper_functions.R \
#            nmt_reach_covariates.tsv surfaces.tsv constants.tsv
args <- commandArgs(trailingOnly = TRUE)
src <- readLines(args[1])
defs <- c("SCALE.vector", "Monomolecular", "Inv.Monomolecular",
          "Inv.Rev.Monomolecular", "Rev.Monomolecular", "Ricker",
          "Inv.Ricker", "Inv.Rev.Ricker", "Rev.Ricker")
env <- new.env()
for (e in parse(text = src)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% defs) eval(e, env)
}
reaches <- read.delim(args[2])

# Parameter sets: variable, transformation, shape, max (ResistanceGA turns
# Ricker-family transformations with shape > 5 into distance, so those
# shapes are at most 5)
sets <- list(
  data.frame(variable = "dor_pc_pva", transformation = "Inv.Monomolecular",
             shape = 1, max = 120),
  data.frame(variable = c("dor_pc_pva", "sgr_dk_rav"),
             transformation = c("Monomolecular", "Ricker"),
             shape = c(3, 2), max = c(50, 200)),
  data.frame(variable = c("dor_pc_pva", "sgr_dk_rav", "ele_mt_cav"),
             transformation = c("Inv.Ricker", "Rev.Monomolecular",
                                "Inv.Rev.Monomolecular"),
             shape = c(4, 8, 2), max = c(10, 700, 3)),
  data.frame(variable = c("sgr_dk_rav", "ele_mt_cav"),
             transformation = c("Rev.Ricker", "Inv.Rev.Ricker"),
             shape = c(5, 3), max = c(1000, 0.5))
)

surfaces <- list()
constants <- list()
for (i in seq_along(sets)) {
  s <- sets[[i]]
  total <- 0
  for (j in seq_len(nrow(s))) {
    x <- env$SCALE.vector(reaches[[s$variable[j]]], 0, 10)
    f <- env[[s$transformation[j]]]
    r <- f(x, parm = c(1, s$shape[j], s$max[j]))
    a <- max(f(x, parm = c(1, s$shape[j], 1))) - 1
    total <- total + r
    constants[[length(constants) + 1]] <- cbind(set = i, s[j, ], A = a)
  }
  surfaces[[i]] <- data.frame(set = i, HYRIV_ID = reaches$HYRIV_ID,
                              value = total / min(total))
}
write.table(do.call(rbind, surfaces), args[3], sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, constants), args[4], sep = "\t",
            row.names = FALSE, quote = FALSE)
