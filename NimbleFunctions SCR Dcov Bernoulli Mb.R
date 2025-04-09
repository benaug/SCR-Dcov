dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0), InSS = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(InSS==1){
      logProb <- log(pi.cell)
    }else{
      logProb <- -Inf
    }
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0), InSS = double(0)) {
    returnType(double(0))
    return(0)
  }
)

GetKern <- nimbleFunction(
  run = function(s = double(1), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)
GetDetectionProb <- nimbleFunction(
  run = function(kern=double(1),p0=double(0),J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      ans <- p0*kern
      return(ans)
    }
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K1D = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){#skip calculation if z=0
      return(0)
    }else{
      logProb <- sum(dbinom(x, size = K1D, p = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0), pd = double(1), K1D = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#Required custom update for number of calls
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("N","y.p","y.c","kern","pd.p","pd.c")) #nodes to copy back to mvSaved
    inds.detected <- control$inds.detected
    z.ups <- control$z.ups
    J <- control$J
    M <- control$M
    #nodes used for update, calcNodes + z nodes
    y.p.nodes <- model$expandNodeNames("y.p")
    y.c.nodes <- model$expandNodeNames("y.c")
    kern.nodes <- model$expandNodeNames("kern")
    pd.p.nodes <- model$expandNodeNames("pd.p")
    pd.c.nodes <- model$expandNodeNames("pd.c")
    N.node <- model$expandNodeNames("N")
    z.nodes <- model$expandNodeNames("z")
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y.p <- model$getLogProb(y.p.nodes[pick])
          lp.initial.y.c <- model$getLogProb(y.c.nodes[pick])

          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0

          #turn kern/pd off
          model$calculate(kern.nodes[pick])
          model$calculate(pd.p.nodes[pick])
          model$calculate(pd.c.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y.p <- model$calculate(y.p.nodes[pick]) #will always be 0
          lp.proposed.y.c <- model$calculate(y.c.nodes[pick]) #will always be 0

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.p + lp.proposed.y.c) - (lp.initial.N + lp.initial.y.p + lp.initial.y.c)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["kern",1][pick,] <<- model[["kern"]][pick,]
            mvSaved["pd.p",1][pick,] <<- model[["pd.p"]][pick,]
            mvSaved["pd.c",1][pick,] <<- model[["pd.c"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["kern"]][pick,] <<- mvSaved["kern",1][pick,]
            model[["pd.p"]][pick,] <<- mvSaved["pd.p",1][pick,]
            model[["pd.c"]][pick,] <<- mvSaved["pd.c",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.p.nodes[pick])
            model$calculate(y.c.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y.p <- model$getLogProb(y.p.nodes[pick]) #will always be 0
          lp.initial.y.c <- model$getLogProb(y.c.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn kern/pd on
          model$calculate(kern.nodes[pick])
          model$calculate(pd.p.nodes[pick])
          model$calculate(pd.c.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y.p <- model$calculate(y.p.nodes[pick])
          lp.proposed.y.c <- model$calculate(y.c.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.p + lp.proposed.y.c) - (lp.initial.N + lp.initial.y.p + lp.initial.y.c)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["kern",1][pick,] <<- model[["kern"]][pick,]
            mvSaved["pd.p",1][pick,] <<- model[["pd.p"]][pick,]
            mvSaved["pd.c",1][pick,] <<- model[["pd.c"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["kern"]][pick,] <<- mvSaved["kern",1][pick,]
            model[["pd.p"]][pick,] <<- mvSaved["pd.p",1][pick,]
            model[["pd.c"]][pick,] <<- mvSaved["pd.c",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.p.nodes[pick])
            model$calculate(y.c.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # copy(from = model, to = mvSaved, row = 1, nodes = z.nodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)