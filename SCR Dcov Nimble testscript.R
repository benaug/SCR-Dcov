library(nimble)
library(coda)
source("sim.SCR.Dcov.R")
source("NimbleModel SCR Dcov Poisson.R")
source("NimbleFunctions SCR Dcov Poisson.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
lam0=0.25
sigma=0.50
K=5
buff=3 #state space buffer. Should be at least 3 sigma.
X<- as.matrix(expand.grid(3:11,3:11))
obstype="poisson"

### Habitat Covariate stuff###
xlim=range(X[,1])+c(-buff,buff)
ylim=range(X[,2])+c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
x.shift=xlim[1]
y.shift=ylim[1]
xlim=xlim-x.shift
ylim=ylim-y.shift
X[,1]=X[,1]-x.shift
X[,2]=X[,2]-y.shift

res=0.25 #habitat grid resolution, length of 1 cell side
cellArea=res^2
x.vals=seq(xlim[1],xlim[2],by=res)
y.vals=seq(ylim[1],ylim[2],by=res)
dSS=as.matrix(expand.grid(x.vals,y.vals))
cells=matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells=nrow(dSS)
n.cells.x=length(x.vals)
n.cells.y=length(y.vals)
#create a density covariate
D.cov=rep(NA,n.cells)
for(c in 1:n.cells){
  D.cov[c]= 6*dSS[c,1] - 0.5*dSS[c,1]^2 + 6*dSS[c,2] - 0.5*dSS[c,2]^2
}
D.cov=as.numeric(scale(D.cov))

plot(dSS,pch=".")
points(X,pch=4,cex=0.75)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Covariate Value")
points(X,pch=4,cex=0.75,col="lightblue")


D.beta0=-2
D.beta1=2
#what is implied expected N in state space?
lambda.cell=exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density")
points(X,pch=4,cex=0.75)

#Simulate some data
data=sim.SCR.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,
                  lam0=lam0,sigma=sigma,K=K,obstype=obstype,
                  X=X,xlim=xlim,ylim=ylim,res)

points(data$s,pch=16)

#Data augmentation level
M=125

#trap operation matrix
J=nrow(X)
K1D=rep(K,J)

#Augment and initialize
n=data$n
y2D=matrix(0,M,J)
y2D[1:n,]=apply(data$y,c(1,2),sum)
xlim=data$xlim
ylim=data$ylim
s.init<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx=which(rowSums(y2D)>0) #switch for those actually caught
for(i in idx){
  trps<- matrix(X[y2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,]<- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,]<- trps
  }
}

z.init=1*(rowSums(y2D)>0)
z.data=rep(NA,M)
z.data[1:n]=1

#inits for nimble
Niminits <- list(z=z.init,s=s.init,lam0=lam0,sigma=sigma,N=sum(z.init>0),D.beta0=0,D.beta1=0)

#constants for Nimble
constants<-list(M=M,J=J,K1D=K1D,xlim=xlim,ylim=ylim,
                D.cov=data$D.cov,cells=cells,cellArea=data$cellArea,n.cells=data$n.cells,
                n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,res=data$res)

#supply data to nimble
dummy.data=rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata<-list(y=y2D,z=z.data,X=X,dummy.data=dummy.data)

# set parameters to monitor
parameters<-c('N','lambda.N','lam0','sigma','D.beta0',"D.beta1")
parameters2 <- c("lambda.cell","s.cell")
nt=1 #thinning rate
nt2=5

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2)


###*required* sampler replacement
z.ups=round(M*0.25) # how many z proposals per iteration???
conf$removeSampler("N")
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(inds.detected=1:n,z.ups=z.ups,J=J,M=M),
                silent = TRUE)


#Nimble-assigned s sampler is fine, but there are two more options here:
#1) update x and y locations jointly using RW_block
#2) update x and y locations jointly, MH for z=1 inds, propose from prior for z=0 inds
#not sure which is faster/yields better mixing.
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=FALSE),silent = TRUE)
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  # type = 'sSampler',control=list(i=i,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
  #                                xlim=data$xlim,ylim=data$ylim,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# Run the model.
start.time2<-Sys.time()
Cmcmc$run(10000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
burnin=250
plot(mcmc(mvSamples[burnin:nrow(mvSamples),]))

#truth
data$N
data$lambda

mvSamples2 = as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx=grep("lambda.cell",colnames(mvSamples2))
burnin2=100


#compare expected D plot to truth
#posterior means
lambda.cell.ests=colMeans(mvSamples2[burnin2:nrow(mvSamples2),lambda.cell.idx])

par(mfrow=c(1,1),ask=FALSE)
zlim=range(c(lambda.cell,lambda.cell.ests)) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)

#cell ests vs. truth
plot(lambda.cell.ests~lambda.cell,pch=16)
abline(0,1,col="darkred",lwd=2)
