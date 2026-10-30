# samc's own model built from rasters (resistance, absorption; fidelity 0;
# conductance 1/mean(x) between 4-neighbours) for the small river of
# test_lattice.py drawn one cell wide on a raster: writes its transient
# matrix Q (with each location's row and column) and the conditional
# passage times among the sites (row = origin, column = destination).
# Usage: Rscript samc_raster.R <samc library> <reference dir>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages({library(Matrix); library(samc); library(terra)})
d <- args[2]
read <- function(f) {
  as.matrix(read.delim(file.path(d, f), header = FALSE))
}
res <- read("samc_raster_resistance.tsv")
abs <- read("samc_raster_absorption.tsv")
fid <- abs * 0
sites <- read.delim(file.path(d, "samc_raster_sites.tsv"))
invisible(capture.output(
  s <- samc(res, abs, fid,
            model = list(fun = "1/mean(x)", dir = 4, sym = TRUE))))
# Cell centres of samc's raster of a matrix (extent 0.5 to n + 0.5, as in
# its R/rasterize.R): x = col, y = nrow - row + 1
xy <- function(r, c) data.frame(x = c, y = nrow(res) - r + 1)
cells <- which(!is.na(res), arr.ind = TRUE)
loc <- locate(s, xy(cells[, 1], cells[, 2]))
Q <- as.matrix(s$q_matrix)
ord <- order(loc)
out <- data.frame(row = cells[ord, 1], col = cells[ord, 2], Q[sort(loc), ])
write.table(out, file.path(d, "samc_raster_Q.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE)
sloc <- locate(s, xy(sites$row, sites$col))
T <- sapply(sloc, function(j) as.numeric(cond_passage(s, dest = j))[sloc])
write.table(T, file.path(d, "samc_raster_cfpt.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE)
cat("samc", as.character(packageVersion("samc")), "\n")
