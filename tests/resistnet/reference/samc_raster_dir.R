# samc's model with a directional transition function, built from rasters
# for the small river of test_lattice.py drawn one cell wide on a raster:
# the transition function of the directional potential of the two cells,
# exp(-(x_to - x_from)) (samc passes the values of the cell moved from,
# then to; sym = FALSE), absorption at the sites, fidelity 0,
# 4-neighbours. Writes its transient matrix Q (with each location's row
# and column) and the conditional passage times among the sites (row =
# origin, column = destination).
# Usage: Rscript samc_raster_dir.R <samc library> <reference dir>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages({library(Matrix); library(samc); library(terra)})
d <- args[2]
read <- function(f) {
  as.matrix(read.delim(file.path(d, f), header = FALSE))
}
pot <- read("samc_raster_potential.tsv")
abs <- read("samc_raster_absorption.tsv")
fid <- abs * 0
sites <- read.delim(file.path(d, "samc_raster_sites.tsv"))
invisible(capture.output(
  s <- samc(pot, abs, fid,
            model = list(fun = function(x) exp(-(x[2] - x[1])), dir = 4,
                         sym = FALSE))))
# Cell centres of samc's raster of a matrix (extent 0.5 to n + 0.5, as in
# its R/rasterize.R): x = col, y = nrow - row + 1
xy <- function(r, c) data.frame(x = c, y = nrow(pot) - r + 1)
cells <- which(!is.na(pot), arr.ind = TRUE)
loc <- locate(s, xy(cells[, 1], cells[, 2]))
Q <- as.matrix(s$q_matrix)
ord <- order(loc)
out <- data.frame(row = cells[ord, 1], col = cells[ord, 2], Q[sort(loc), ])
write.table(out, file.path(d, "samc_raster_dir_Q.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE)
sloc <- locate(s, xy(sites$row, sites$col))
T <- sapply(sloc, function(j) as.numeric(cond_passage(s, dest = j))[sloc])
write.table(T, file.path(d, "samc_raster_dir_cfpt.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE)
cat("samc", as.character(packageVersion("samc")), "\n")
