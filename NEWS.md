geonet version 0.3.0 (Github release)

- algorithm options can now be supplied to intensity_pspline via the
  "control" argument
- as_gn now also takes the "units" attribute from a linnet object  
- added missing summary and print methods
- fixed bug in "internal" occurring when a linear covariates has length 1
- internal covariates are recognized automatically and do not need to be 
  assigned in the formula

geonet version 0.2.0 (GitHub release)

- added internal covariate information to the montgomery network
- fixed bug in network_penalty occurring when r = 1
- model summary has been changed
- fit_poisson_model is only doing on Fisher scoring iteration within the inner 
  loop from the second iteration
- effective degrees of freedom is computed for smooth terms
- allow for general internal covariates added to the lins attribute of the
  network
- internal covariate "dist2V" can be added to the linear predictor 
- stopping criterion now depends on relative difference of theta and not of rho
- verbose argument added to intensity_pspline which allows to track the 
  progress of the fitting algorithm

geonet version 0.1.1 (CRAN and GitHub release)

- Added reference in the description file.
- Added missing return values in the documentation.
- Removed set.seed() call.

geonet version 0.1.0 (no release)

- Initial submission
