# Lundgren & Ralph's (2019) own HCalc (gene.flow.inference, R/mcmc.fns.R)
# for the small river graph of test_coalescent.py: the mean coalescence
# times (symmetrised, as their H_and_lllh), written to lundgren_ralph_H.tsv.
# Their package has no licence, so its code is not copied here; pass the
# path of its R/ directory (github.com/elundgre/gene-flow-inference).
# Usage: Rscript lundgren_ralph.R <gene.flow.inference/R> <reference dir>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(Matrix))
source(file.path(args[1], "mcmc.fns.R"))
d <- args[2]
edges <- read.delim(file.path(d, "lundgren_ralph_edges.tsv"))
gam <- scan(file.path(d, "lundgren_ralph_gamma.txt"), quiet = TRUE)
n <- length(gam)
G_adj <- sparseMatrix(i = edges$x + 1, j = edges$y + 1, x = 1,
                      dims = c(n, n))
# G_adj@x is column by column: put the rates in that order
g <- edges$rate[order(edges$y, edges$x)]
H <- HCalc(g, gam, G_adj, n, const_coal = FALSE)
H <- (H + t(H)) / 2
write.table(as.matrix(H), file.path(d, "lundgren_ralph_H.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE)
