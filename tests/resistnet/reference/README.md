# Reference values for the NMT example data

Computed with independent implementations, so that tests check ResistNet
against them rather than against its own earlier output.

- `samc_raster_*.tsv`: samc's own SAMC (samc 3.2.1), built by samc from
  rasters: a small river of five reaches (`test_lattice.py`) drawn one
  cell wide on a 12 x 16 raster, with each reach's resistance per km and a
  confluence cell taking the mean of the reaches meeting there
  (`samc_raster_resistance.tsv`), absorption at five sites
  (`samc_raster_absorption.tsv`, `samc_raster_sites.tsv`), 4-neighbours and
  conductance 1/mean(x). `samc_raster.R` writes samc's transient matrix
  (`samc_raster_Q.tsv`, each location's row and column first) and the
  conditional passage times among the sites (`samc_raster_cfpt.tsv`).
  `python test_lattice.py` writes the input rasters. The lattice
  (`resistnet/lattice.py`) must give the same matrix and passage times.
- `nmt_samc_lattice_*.tsv`: conditional first passage times among the 14
  sites (row = origin, column = destination), from `samc::cond_passage`
  (samc 3.2.1, `samc(data = P)`), for three SAMC models on the lattice of
  `nmt.gpkg` (see `test_samc_network.py::R_SAMC_MODELS`):
  `nmt_samc_lattice.py` exports their transition matrices and
  `nmt_samc_cond_passage.R` runs samc.
- `nmt_mlpe_Y.tsv`: the NMT genetic distances (linearised FST,
  `nmt.popGenDistMat.txt`), in the network's site order. `nmt_mlpe_X_*.tsv`: model distances among the
  sites: river distance in km (`streamkm`), the standard model with
  degree of regulation (`dor_pc_pva`, inverse monomolecular, shape 1;
  `dor_inv_mono`), and the symmetrised passage times of the `sgr_flow`
  SAMC model on the lattice (`samc_sgr_flow`, from
  `nmt_samc_lattice.py`).
- `nmt_mlpe_lme4.tsv`: ResistanceGA's `MLPE.lmm` fitted to Y and each X
  (lme4, ML; AIC counts the four MLPE parameters only) and the marginal
  R2 from `MuMIn::r.squaredGLMM` (`nmt_mlpe_lme4.R` made it).
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
