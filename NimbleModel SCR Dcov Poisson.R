NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  #Dcov priors
  D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #detection priors
  lam0 ~ dunif(0,10)
  sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  #Density model
  for(c in 1:n.cells) {
    lambda.cell[c] <- exp(D.beta0 + D.beta1*D.cov[c])*cellArea #expected N in cell c
    pi.cell[c] <- lambda.cell[c] / lambda.N
  }
  lambda.N <- sum(lambda.cell[1:n.cells])
  N ~ dpois(lambda.N)
  
  for(i in 1:M) {
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~  dunif(xlim[1],xlim[2])
    s[i,2] ~  dunif(ylim[1],ylim[2])
    #get cell s i lives in
    s.cell[i] <- getCell(s=s[i,1:2],res=res,cells=cells[1:n.cells.x,1:n.cells.y])
    dummy.data[i] ~ dCell(s.cell[i],pi.cell[1:n.cells]) #get categorical density likelihood for this cell
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
})