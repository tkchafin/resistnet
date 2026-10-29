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
- `resistancega_transformations.tsv`: ResistanceGA's eight transformation
  functions (`Ricker`, `Rev.Ricker`, `Inv.Ricker`, `Inv.Rev.Ricker`,
  `Monomolecular`, `Rev.Monomolecular`, `Inv.Monomolecular`,
  `Inv.Rev.Monomolecular`, vector versions from its
  `R/internal_helper_functions.R`) applied in R to a test covariate
  rescaled to 0-10 with its `SCALE.vector`, for shapes 0.5-50, then
  rescaled to 0-1 (`resistancega_transformations.R` made it). Note that
  ResistanceGA itself turns any Ricker-family transformation with shape
  > 5 into plain distance; ResistNet does not.
- `resistancega_surfaces.tsv`: ResistanceGA's composite resistance
  surface (as in its `R/Combine_Surfaces.R`: each continuous surface
  rescaled to 0-10 with `SCALE.vector`, transformed with parm =
  c(eq, shape, max), summed, and divided by its minimum) for the 127 NMT
  reaches of `nmt.network`, for four sets of one to three surfaces (six
  of the eight transformations; maxima from 0.5 to 1000), computed in R
  with its own functions (`resistancega_surfaces.R`, from
  `nmt_reach_covariates.tsv`, the reaches' raw values).
  `resistancega_surface_constants.tsv` gives each surface's parameters
  and its range at max = 1 (A), which map them onto ResistNet's
  resistance per km (weights max_i x A_i, and K - 1 = sum(max_i x A_i) /
  number of surfaces).
