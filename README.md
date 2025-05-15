# SCR-Dcov
Nimble MCMC samplers for Spatial Capture Recapture with density covariates.

Two unique features:

1) Use of continuous state space with discrete grid for habitat covariates. This should be faster than typical discrete state space approaches.
2) These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation


Some additional notes:
1) Care should be taken with state spaces with multiple habitat patches surrounded by non-habitat, or even very low D habitat. The concern is that activity centers will be able to move within a habitat patch during sampling, but not across patches, once the random walk proposal is tuned. Using the provided "sSampler" bypasses this problem for uncaptured and augmented individuals by proposing from the combined spatial prior (categorical for discrete state space cell, then uniform within cell) when z=0. This allows activity centers to move across patches when their z_i=0 and then be turned on/off in the new patches. Therefore, this sampler for activity centers should be preferred over a regular random walk update.

The remaining concern is whether any *captured* individuals are stuck in one habitat patch when they should be moving between them. This seems like a rare occurrence in practice, but I cannot rule out that it may happen with this approach and introduce some degree of bias and/or uncertainty underestimation. If this happens, it is probably best to treat the state space as completely discrete as is typically done.