library(nimble)
library(coda)
source("sim.SCR.Dcov.R")
source("NimbleModel SCR Dcov Bernoulli.R")
source("NimbleFunctions SCR Dcov Bernoulli.R")
source("sSampler.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

#simulate some data
p0 <- 0.25 #baseline detection probability
sigma <- 0.50 #detection spatial scale
K <- 5 #number of occasions
buff <- 2 #state space buffer. Should be at least 3 sigma (generally).
X <- as.matrix(expand.grid(3:11,3:11)) #trapping array
obstype <- "bernoulli" #observation model

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1])+c(-buff,buff)
ylim <- range(X[,2])+c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#create a density covariate - one for each session
library(geoR)
#need a simulated landscape with individuals living around traps to be captured
set.seed(13225)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(100,100),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale

par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4,cex=0.75,col="darkred",lwd=2)

#Additionally, maybe we want to exclude "non-habitat" or limit the state space extent
#let's use a 3sigma buffer
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(0,length(D.cov))
dists <- e2dist(X,dSS.tmp)
min.dists <- apply(dists,2,min)
InSS[min.dists<(3*sigma)] <- 1
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="Habitat",col=cols1)
points(X,pch=4,col="darkred",lwd=2)

#Density covariates
D.beta0 <- -1
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density")
points(X,pch=4,cex=0.75)

#Simulate some data
set.seed(33923) #seed set above for D.cov. Change seed here to simulate different data sets
data <- sim.SCR.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                  p0=p0,sigma=sigma,K=K,obstype=obstype,
                  X=X,xlim=xlim,ylim=ylim,res=res)

points(data$s,pch=16)

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
                       n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
                       x.vals=data$x.vals,y.vals=data$y.vals)

#Data augmentation level
M <- 150

#trap operation matrix
J <- nrow(X)
K1D <- rep(K,J)

#Augment and initialize
n <- data$n
y2D <- matrix(0,M,J)
y2D[1:n,] <- apply(data$y,c(1,2),sum)
xlim <- data$xlim
ylim <- data$ylim
s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx <- which(rowSums(y2D)>0) #switch for those actually caught
for(i in idx){
  trps <- matrix(X[y2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,] <- trps
  }
}

#If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
e2dist  <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
getCell  <-  function(s,res,cells){
  cells[trunc(s[1]/data$res)+1,trunc(s[2]/data$res)+1]
}
alldists <- e2dist(s.init,data$dSS)
alldists[,data$InSS==0] <- Inf
for(i in 1:M){
  this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
  if(data$InSS[this.cell]==0){
    cands <- alldists[i,]
    new.cell <- which(alldists[i,]==min(alldists[i,]))
    s.init[i,] <- data$dSS[new.cell,]
  }
}

#plot to make sure initialized activity centers are in habitat
image(data$x.vals,data$y.vals,matrix(data$InSS,data$n.cells.x,data$n.cells.y))
points(s.init,pch=16)


z.init <- 1*(rowSums(y2D)>0)
z.data <- rep(NA,M)
z.data[1:n] <- 1

#inits for nimble - MUST use z init and N init for data augmentation scheme to work. should use s.init, too.
Niminits <- list(z=z.init,N=sum(z.init>0),s=s.init,p0=runif(1,0.2,0.8),sigma=runif(1,0.2,0.8),
                 D0=sum(z.init)/(sum(data$InSS)*data$res^2),D.beta1=0)

#constants for Nimble
#here, you probably want to center your D.cov. The one I simulated for this testscript is already centered.
# D.cov.use <- data$D.cov - mean(data$D.cov) #plug this into constants$D.cov if centering
constants <- list(M=M,J=J,K1D=K1D,xlim=xlim,ylim=ylim,
                D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                res=data$res)

#supply data to nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y=y2D,z=z.data,X=X,dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)

# set parameters to monitor
parameters<-c('N','lambda.N','p0','sigma','D0',"D.beta1")
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 5 #thinning rate for paremeters2

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c("p0","sigma")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2,thin2=nt2,
                      nodes=config.nodes,useConjugacy=FALSE)

#add N/z sampler
z.ups <- round(M*0.25) # how many z proposals per iteration??? 25% of M generally seems good, but no idea what is optimal
conf$removeSampler("N")
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(inds.detected=1:n,z.ups=z.ups,J=J,M=M),
                silent = TRUE)

#add s sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                 xlim=data$xlim,ylim=data$ylim),silent = TRUE)
}

#add blocked sampler for Dcovs
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
#add block sampler for p0 and sigma. AF_slice slow here, using RW_block
conf$addSampler(target = c("p0","sigma"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[burnin:nrow(mvSamples),]))

#truth
data$N
data$lambda

mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10

#compare expected D plot to truth
#image will show 
#posterior means
lambda.cell.post <- cellArea*mvSamples2[burnin2:nrow(mvSamples2),D0.idx]*mvSamples2[burnin2:nrow(mvSamples2),lambda.cell.idx]
lambda.cell.ests <- colMeans(lambda.cell.post)
#remove non-habitat
lambda.cell.ests[InSS==0] <- NA
lambda.cell[InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)