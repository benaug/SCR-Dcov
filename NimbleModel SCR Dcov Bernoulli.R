NimModel <- nimbleCode({
  #Density parameter priors
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #detection priors
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  #Density model
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space
  
  for(i in 1:M) {
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~  dunif(xlim[1],xlim[2])
    s[i,2] ~  dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InSS=InSS[s.cell[i]])
    #Observation model, skipping z_i=0 calculations
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y[i,1:J] ~ dBinomialVector(pd=pd[i,1:J],K1D[1:J],z=z[i]) #vectorized obs mod
  }
})
#custom Metropolis-Hastings update for N/z
