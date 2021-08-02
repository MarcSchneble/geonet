geonet version 0.6.0 (Github release)

- revised plotting functions (caution: new arguments defining line widths, etc.)
- scale of the dist2V variable is now supported 
- added argument "other arguments" to intensity_pspline for compatibility with older versions
- revised summary for gnppfit objects

geonet version 0.5.0 (Github release)

- new function network_ISE for the computation of integrated squared errors
- new function network_integral for the computation of an integral of the intensity
- geonet supplies wrappers for density.lpp, i.e. the intensity on a geometric network can be
  estimated by employing kernel based methods
- fixed bug which occurs if a smooth covariate is specified in formula after a linear covariate

geonet version 0.4.0 (Github release)

- summaries and related print methods now also print information on the
  linear network representation
- default penalty is now of order r = 2
- the global knot distance delta can now be supplied as a quantile of the curve
  lenghts of the network and the global bin width h as a fraction of the global
  knot distance
- fixed bug in as_gnpp.lpp when linear point pattern is not marked
- fixed bug in runif_gn
- added function rgnpp which allows to simulate from a fitted intensity

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
