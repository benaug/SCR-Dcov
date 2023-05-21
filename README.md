# SCR-Dcov
Nimble MCMC samplers for Spatial Capture Recapture with density covariates.

Two unique features:

1) Use of continuous state space with discrete grid for habitat covariates. This should be faster than typical discrete state space approaches.
2) Alternative data augmentation scheme, similar to reversible jump MCMC, that allows realized abundance to be Poisson.
