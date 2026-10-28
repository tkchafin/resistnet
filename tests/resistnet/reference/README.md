# Reference values for the NMT example data

Computed with independent implementations, so that tests check ResistNet
against them rather than against its own earlier output.

- `nmt_samc_*.tsv`: conditional first passage times among the 14 sites
  (row = origin, column = destination), from `samc::cond_passage` in R
  (samc 3.2.1), for three SAMC models on `nmt.gpkg` (see
  `test_samc_network.py::R_SAMC_MODELS`). They also agree with a 60-digit
  decimal solution (dense Gaussian elimination) to 2e-13.
- `nmt_mlpe_Y.tsv`: the NMT genetic distances (pairwise FST), in the
  network's site order. `nmt_mlpe_X_*.tsv`: model distances among the
  sites: river distance in km (`streamkm`), the standard model with
  degree of regulation (`dor_pc_pva`, inverse monomolecular, shape 1;
  `dor_inv_mono`), and the symmetrised passage times of the
  `ele_dir_flow` SAMC model.
- `nmt_mlpe_lme4.tsv`: ResistanceGA's `MLPE.lmm` fitted to Y and each X
  (lme4, ML; AIC counts the four MLPE parameters only) and the marginal
  R2 from `MuMIn::r.squaredGLMM`.
