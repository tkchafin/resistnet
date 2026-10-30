# Lundgren & Ralph's (2019) own sampler (gene.flow.inference, R/findG.MH.R)
# on the divergences of lundgren_ralph.R (same graph and coalescence rates
# per deme, default error variance, adaptive proposals), from the true
# rates, with its random numbers replaced by those of
# lundgren_ralph_mh_normals.txt (proposals) and lundgren_ralph_mh_uniforms.txt
# (acceptances, as u < p), so that resistnet.coalescent.mh() can be fed the
# same ones and must take the same steps (test_coalescent.py). Writes the
# chain (rates in the order of lundgren_ralph_edges.tsv), log posterior,
# rejections and proposal standard deviations at every iteration. Their
# package has no licence, so its code is not copied here; pass the path of
# its R/ directory (github.com/elundgre/gene-flow-inference).
# Usage: Rscript lundgren_ralph_mh.R <gene.flow.inference/R> <reference dir>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(Matrix))
source(file.path(args[1], "mcmc.fns.R"))
source(file.path(args[1], "findG.MH.R"))
d <- args[2]
edges <- read.delim(file.path(d, "lundgren_ralph_edges.tsv"))
gam <- scan(file.path(d, "lundgren_ralph_gamma.txt"), quiet = TRUE)
D <- as.matrix(read.delim(file.path(d, "lundgren_ralph_inputs.tsv"),
                          header = FALSE))
z <- scan(file.path(d, "lundgren_ralph_mh_normals.txt"), quiet = TRUE)
u <- scan(file.path(d, "lundgren_ralph_mh_uniforms.txt"), quiet = TRUE)
# findG.MH is defined here, so it finds these before the stats package's
iz <- 0
iu <- 0
rnorm <- function(n, mean = 0, sd = 1) {
  x <- z[iz + seq_len(n)]
  iz <<- iz + n
  mean + sd * x
}
rbinom <- function(n, size, prob) {
  x <- u[iu + seq_len(n)]
  iu <<- iu + n
  as.integer(x < prob)
}
n <- length(gam)
iter <- length(z) / (nrow(edges) + n)
G_adj <- sparseMatrix(i = edges$x + 1, j = edges$y + 1, x = 1,
                      dims = c(n, n))
# G_adj@x is column by column: rates in that order, and back
o <- order(edges$y, edges$x)
ans <- findG.MH(D, G_adj, const_coal = FALSE, iter = iter,
                fixed_start = TRUE, g_init = edges$rate[o], gam_init = gam)
g <- ans$g[, order(o), drop = FALSE]
colnames(g) <- paste0("g_", edges$x, "_", edges$y)
gm <- ans$gam
colnames(gm) <- paste0("gam_", seq_len(n) - 1)
out <- data.frame(g, gm, lp = ans$lllh, reject_g = ans$reject_g,
                  reject_gam = ans$reject_gam, sd_g = ans$sdG,
                  sd_gam = ans$sdgam)
write.table(format(out, digits = 17), file.path(d, "lundgren_ralph_mh.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
