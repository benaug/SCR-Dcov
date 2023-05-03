# SCR-Dcov
Nimble MCMC samplers for Spatial Capture Recapture with density covariates.

Two unique features:

1) Use of continuous state space with discrete grid for habitat covariates. This should be faster than typical discrete state space approaches.
2) Alternative data augmentation scheme, similar to reversible jump MCMC, that allows realized abundance to be Poisson.


Disclaimer: This has not been tested, yet, but looks plausibly correct from a few runs. I will test this sampler via simulation and then and remove the disclaimer.