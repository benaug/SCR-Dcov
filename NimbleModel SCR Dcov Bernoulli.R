NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  #Dcov priors
  D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #detection priors
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  #Density model
  for(c in 1:n.cells) {
    lambda.cell[c] <- InHabitat[c]*exp(D.beta0 + D.beta1*D.cov[c])*cellArea #expected N in cell c
    pi.cell[c] <- lambda.cell[c] / lambda.N #expected proportion of total N in cell c
  }
  lambda.N <- sum(lambda.cell[1:n.cells]) #expected N in state space
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
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InHabitat=InHabitat[s.cell[i]])
    #Observation model, skipping z_i=0 calculations
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y[i,1:J] ~ dBinomialVector(pd=pd[i,1:J],K1D[1:J],z=z[i]) #vectorized obs mod
  }
})
#custom Metropolis-Hastings update for N/z
