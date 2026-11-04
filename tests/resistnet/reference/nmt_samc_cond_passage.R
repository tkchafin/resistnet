# samc::cond_passage among the NMT sites for each lattice SAMC model
# exported by nmt_samc_lattice.py (samc(data = P), P with an absorbing
# state): writes nmt_samc_lattice_<model>.tsv (row = origin, column =
# destination) in this directory.
# Usage: Rscript nmt_samc_cond_passage.R <samc library> <scratch dir>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages({library(Matrix); library(samc)})
for (f in list.files(args[2], pattern = "\\.P\\.mtx$", full.names = TRUE)) {
  name <- sub("\\.P\\.mtx$", "", basename(f))
  P <- as(as(readMM(f), "CsparseMatrix"), "generalMatrix")
  invisible(capture.output(s <- samc(data = P)))
  dests <- scan(file.path(args[2], paste0(name, ".dests.txt")),
                quiet = TRUE)
  T <- sapply(dests, function(d) as.numeric(cond_passage(s, dest = d))[dests])
  write.table(T, paste0("nmt_samc_lattice_", name, ".tsv"), sep = "\t",
              row.names = FALSE, col.names = FALSE)
}
cat("samc", as.character(packageVersion("samc")), "\n")
