library(nimble)
library(coda)
source("sim.SCR.Dcov.R")
source("NimbleModel SCR Dcov Bernoulli.R")
source("NimbleFunctions SCR Dcov Bernoulli.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
p0 <- 0.25 #baseline detection probability
sigma <- 0.50 #detection spatial scale
K <- 5 #number of occasions
buff <- 3 #state space buffer. Should be at least 3 sigma (generally).
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

res <- 0.20 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1],xlim[2],by=res)
y.vals <- seq(ylim[1],ylim[2],by=res)
dSS <- as.matrix(expand.grid(x.vals,y.vals)) + res/2 #add res/2 to get cell centroids
#remove extra cells created outside xlim and ylim
rem.idx <- which(dSS[,1]>xlim[2]|dSS[,2]>ylim[2])
dSS <- dSS[-rem.idx,]
cells <- matrix(1:nrow(dSS),nrow=length(x.vals)-1,ncol=length(y.vals)-1)
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals) - 1
n.cells.y <- length(y.vals) - 1

#create a density covariate
D.cov <- rep(NA,n.cells)
for(c in 1:n.cells){
  D.cov[c] <-  7*dSS[c,1] - 0.5*dSS[c,1]^2 + 7*dSS[c,2] - 0.5*dSS[c,2]^2
}
D.cov <- as.numeric(scale(D.cov))

#Visualize dSS with grid. dSS should be grid centroids
plot(dSS,pch=".")
abline(v=x.vals,h=y.vals)

image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Covariate Value")
points(X,pch=4,cex=0.75,col="lightblue")


#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InHabitat=rep(1,length(D.cov))
InHabitat[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InHabitat[dSS.tmp[,1]<2&dSS.tmp[,2]>12] <- 0
InHabitat[dSS.tmp[,1]>12&dSS.tmp[,2]<2] <- 0
InHabitat[dSS.tmp[,1]>12&dSS.tmp[,2]>12] <- 0

image(x.vals,y.vals,matrix(InHabitat,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -2
D.beta1 <- 2
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density")
points(X,pch=4,cex=0.75)

#Simulate some data
data <- sim.SCR.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InHabitat=InHabitat,
                  p0=p0,sigma=sigma,K=K,obstype=obstype,
                  X=X,xlim=xlim,ylim=ylim,res=res)

points(data$s,pch=16)

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
s.init<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx <- which(rowSums(y2D)>0) #switch for those actually caught
for(i in idx){
  trps<- matrix(X[y2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,]<- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,]<- trps
  }
}

#If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
e2dist  <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
getCell  <-  function(s,res,cells){
  cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
}
alldists <- e2dist(s.init,data$dSS)
alldists[,data$InHabitat==0] <- Inf
for(i in 1:M){
  this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
  if(data$InHabitat[this.cell]==0){
    cands <- alldists[i,]
    new.cell <- which(alldists[i,]==min(alldists[i,]))
    s.init[i,] <- data$dSS[new.cell,]
  }
}

#plot to make sure initialized activity centers are in habitat
image(data$x.vals,data$y.vals,matrix(data$InHabitat,data$n.cells.x,data$n.cells.y))
points(s.init,pch=16)


z.init <- 1*(rowSums(y2D)>0)
z.data <- rep(NA,M)
z.data[1:n] <- 1

#inits for nimble - MUST use z init and N init for data augmentation scheme to work. should use s.init, too.
Niminits <- list(z=z.init,N=sum(z.init>0),s=s.init,
                 p0=p0,sigma=sigma,D.beta0=0,D.beta1=0)

#constants for Nimble
constants<-list(M=M,J=J,K1D=K1D,xlim=xlim,ylim=ylim,
                D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                res=data$res)

#supply data to nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata<-list(y=y2D,z=z.data,X=X,dummy.data=dummy.data,cells=cells,InHabitat=data$InHabitat)

# set parameters to monitor
parameters<-c('N','lambda.N','p0','sigma','D.beta0',"D.beta1")
parameters2 <- c("lambda.cell","s.cell")
nt <- 1 #thinning rate for parameters
nt2 <- 5 #thinning rate for paremeters2

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2,thin2=nt2,
                      useConjugacy=FALSE) #configure very slow if checking for conjugacy


###*required* sampler replacement
z.ups <- round(M*0.25) # how many z proposals per iteration??? 25% of M generally seems good, but no idea what is optimal
conf$removeSampler("N")
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(inds.detected=1:n,z.ups=z.ups,J=J,M=M),
                silent = TRUE)

#Can try block RW sampler if D.cov posteriors strongly correlated
# conf$removeSampler(c("D.beta0","D.beta1"))
# conf$addSampler(target = c("D.beta0","D.beta1"),type = 'RW_block',control = list(),silent = TRUE)


#Nimble-assigned s sampler is fine, but there are two more options here:
#1) update x and y locations jointly using RW_block
#2) update x and y locations jointly, MH for z=1 inds, propose from prior for z=0 inds
#not sure which is faster/yields better mixing. Proposing from dSS will be slower with more cells
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  #                 type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=FALSE),silent = TRUE)
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  type = 'sSampler',control=list(i=i,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                 xlim=data$xlim,ylim=data$ylim,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <-as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[burnin:nrow(mvSamples),]))

#truth
data$N
data$lambda

mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
burnin2 <- 10


#compare expected D plot to truth
#image will show 
#posterior means
lambda.cell.ests <- colMeans(mvSamples2[burnin2:nrow(mvSamples2),lambda.cell.idx])
#remove non-habitat
lambda.cell.ests[InHabitat==0] <- NA
lambda.cell[InHabitat==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)

#cell ests vs. truth
plot(lambda.cell.ests~lambda.cell,pch=16) #remove non-habitat
abline(0,1,col="darkred",lwd=2)
