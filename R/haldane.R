#
# 	haldane.R
#
#	$Revision: 1.22 $	$Date: 2008/06/23 05:49:16 $
#
####################################################################
#
#  Calculations for Haldane type models
#

saturated.state <- function(model="D", depth=0, g=air) {
  if(is.character(model))
    model <- pickmodel(model)
  else
    stopifnot(inherits(model, "hm"))
  if(!capable(model, g))
    stop(paste("Model cannot deal with", as.character(g)))
  y <- summary(model)
  Pamb <- 1 + depth/10
  out <- list(N2= rep(Pamb * g$fN2, y$nc), He=rep(Pamb * g$fHe, y$nc))
  out <- out[y$species]
  out <- as.data.frame(out)
  return(out)
}

conform <- function(state, model="D") {
  if(is.character(model))
    model <- pickmodel(model)
  else
  stopifnot(inherits(model, "hm"))
  m <- summary(model)
  mgas <- m$species
  if(is.numeric(state) && is.vector(state))
    s <- data.frame(N2=state)
  else if(is.list(state)) 
    s <- as.data.frame(state)
  else if(is.matrix(state) || is.data.frame(state))
    s <- state
  else
    stop("Unrecognised format for state")
  sgas <- colnames(s)
  if(is.null(sgas)) {
    if(ncol(s) == length(mgas))
      colnames(s) <- mgas
    else
      stop("gases in state do not match gases in model")
  } else {
    if(all(sgas %in% mgas) && all(mgas %in% sgas))
      s <- s[mgas]
    else
      stop("gases in state do not match gases in model")
  }
  if(nrow(s) != m$nc)
    stop("numbers of compartments in state do not match compartments in model")
  return(s)
}

"haldane" <- 
function(d,
	 model=pickmodel("DSAT"),
	 prevstate=NULL,
         progressive=FALSE) 
{
  stopifnot(is.dive(d))

  if(is.character(model))
    model <- pickmodel(model)
  else
    stopifnot(inherits(model, "hm"))
  y <- summary(model)

  if(is.null(prevstate))
    prevstate <- saturated.state(model, 0)
  else
    prevstate <- conform(prevstate, model)
  
  if(!capable(model, d, "M0"))
    stop("The model does not include data for some of the gases in this dive")

  species <- y$species
  
  # vectors: time
  times <- times.dive(d)
  depths <- depths.dive(d)
  durations <- diff(times)
  n <- length(times)

  # data frames: compartment x species
  state <- as.data.frame(prevstate)
  
  HalfT <- lapply(species,
                  function(x, ...) { param(species=x, what="HalfT", ...) },
                  model=model)
  names(HalfT) <- species
  HalfT <- data.frame(HalfT)
  ratecoefs <- log(2)/HalfT

  ncomp <- nrow(state)
  
  # data frame: time x species
  fGas <- with(d$data, list(N2=fN2, He=fHe))
  fGas <- as.data.frame(fGas[species])

  if(progressive) {
    profile <- array(0, dim=c(n, dim(state)))
    profile[1,,] <- as.matrix(state)
  }

  for(i in 1:(n-1)) {
    d0 <- depths[i]
    d1 <- depths[i+1]
    tim <- durations[i]
    if(tim > 0) {
      fgi <- as.numeric(fGas[i,])
      k1 <- fgi * (d0/10 + 1)
      k2 <- fgi * ((d1-d0)/10) / tim
      k1 <- outer(rep(1, ncomp), k1, "*")
      k2 <- outer(rep(1, ncomp), k2, "*")
      v <- exp(- ratecoefs * tim)
      state <- state * v + ((k1 - k2/ratecoefs) * (1 - v) + k2 * tim)
    }
    if(progressive)
      profile[i+1,,] <- as.matrix(state)
  }
  
  rownames(state) <- rownames(model)
  colnames(state) <- species

  if(progressive)
    return(profile)
  else 
    return(state)
}


###################################################################
#
#  Interactive display
#

showstates <- function(d, model="DSAT") {
  stopifnot(is.dive(d))

  # pick model by its name
  if(is.character(model))
    model <- pickmodel(model)
  else if(!inherits(model, "hm"))
    stop("model should be an object of class hm")

  if(!capable(model, d, "M0"))
    stop("The model does not include data for some of the gases in this dive")

  y <- summary(model)

  Mvals <- M0mix(model, d$data$fN2, d$data$fHe)
  
  oldpar <- par(mfrow=c(1,2), ask=FALSE)
  frame()
  plot(d)
  timepoints <- times.dive(d)
  maxtime <- max(timepoints)
  ntissues <- y$nc
  
  # precompute time profile of tissue saturation
  cat("precomputing tissue saturations...")
  halprofile <- haldane(d, model, progressive=TRUE)

  # sum over gas species N2 + He = 'inert gas'
  igprofile <- apply(halprofile, c(1,2), sum) 
                        
  # convert to relative saturation profile
  relprofile <- igprofile/Mvals
  
  # determine (approx) maximum relative saturation over entire dive
  maxrelsat <- max(max(relprofile), 1)

  # precompute oxygen toxicity profile
  oxtoxprofile <- oxtox(d, progressive=TRUE)

  cat("done.\nplease click on the graph\n")

  # event loop ..
  hal <- rep(NA, ntissues)

  while(length(xy <- locator(1)) > 0) {
    tim <- xy$x
    tim <- max(0, tim)
    tim <- min(tim, maxtime)

    ind  <- max(which(timepoints <= tim))
    ind1 <- min(which(timepoints >= tim))
    hal <- halprofile[ind,,] 
    ig  <- igprofile[ind,] 
    rel <- as.numeric(relprofile[ind, ])
    otu <- oxtoxprofile[ind]

    depths <- depths.dive(d)
    if(ind < ind1) {
      # interpolate
      t0 <- timepoints[ind]
      t1 <- timepoints[ind1]
      divebit <- chop.dive(dive.segment(d, ind), t0, tim)
      hal <- haldane(divebit, model, prevstate=hal)
      ig <- apply(hal, 1, sum)
      M0. <- M0mix(model, d$data$fN2[ind], d$data$fHe[ind])
      rel <- as.numeric(ig/M0.)
    }
    pozzie <- barplot(rel,
                      xlab=paste("Tissues (", y$title, ")", sep=""),
                      ylab="Relative saturation",
                      main=paste("Time=", round(tim), "min\n",
                        "Accumulated oxygen toxicity", round(otu, 1)),
                      ylim=range(c(0, 1.1 * maxrelsat)),
                      names.arg=NULL)
    if(ntissues > 1) {
      mtext(side=1, at=pozzie[1], line=1, text="Fast")
      mtext(side=1, at=pozzie[ntissues], line=1, text="Slow")
    }
    abline(h=1, lty=3, col="red")
    plot(d)
    abline(v=tim, lty=2, col="green")
  }
  par(oldpar)
  if(is.matrix(hal))
    colnames(hal) <- y$species
  hal
}


###################################################################
#
#  NDL calculation
#

ndl <- function(depth, g=air, prevstate=NULL, model="DSAT") {
  stopifnot(is.gas(g))
  if(is.character(model))
    model <- pickmodel(model)
  else stopifnot(inherits(model, "hm"))
  if(!capable(model, g))
    stop(paste("Model cannot deal with", as.character(g)))
  if(is.null(prevstate))
    prevstate <- saturated.state(model)
  else
    prevstate <- conform(prevstate, model)
  result <- numeric(n <- length(depth))
  if(n == 0) return(result)
  species <- summary(model)$species
  pressure <- depth/10 + 1
  if(all(species == "N2")) {
    # N2 only
    fN2 <- g$fN2
    M0 <- param(model, "N2", "M0")
    ratecoefs <- log(2)/param(model, "N2", "HalfT")
    init <- prevstate$N2
    for(i in 1:n) {
      ppN2 <- fN2 * pressure[i]
      result[i] <- 
        if(any(init > M0))
          # initial state of diver does not permit a dive
          0
        else if(!any(bite <- (ppN2 > M0)))
        # tissue tension will never exceed M-value in any compartment
          Inf
        else
          min(-log(((ppN2-M0)[bite])/((ppN2-init)[bite]))/ratecoefs[bite])
    }
    return(result)
  } else {
    # N2 and He 
    fN2 <- g$fN2
    fHe <- g$fHe
    rateN2 <- log(2)/param(model, "N2", "HalfT")
    rateHe <- log(2)/param(model, "He", "HalfT")
    rateSLOWER <- pmin(rateN2, rateHe)
    ntissues <- length(rateN2)
    # compute surfacing M-value for gas
    M0 <- M0mix(model, fN2, fHe)
    # initial tissue saturations
    initN2 <- prevstate$N2
    initHe <- prevstate$He
    initIG <- initN2 + initHe
    #
    decay <- function(x, init, rate, asymp) {
      asymp + exp(-rate * x) * (init-asymp)
    }
    objective <- function(x, ppN2, initN2, rateN2, ppHe, initHe, rateHe, M0) {
      decay(x, initN2, rateN2, ppN2) + decay(x, initHe, rateHe, ppHe) - M0
    }
    for(i in 1:n) {
      ppN2 <- fN2 * pressure[i]
      ppHe <- fHe * pressure[i]
      ppIG <- ppN2 + ppHe
      val <- numeric(ntissues)
      for(j in 1:ntissues) {
        if(initIG[j] > M0[j])
          # max tissue tension already exceeded
          val[j] <- 0
        else if(ppIG <= M0[j])
          # tissue tension will never exceed M-value in this compartment
          val[j] <-   Inf
        else {
          endpoint <- -log((ppIG-M0[j])/(ppIG-initIG[j]))/rateSLOWER[j]
          val[j] <- uniroot(objective,
                            c(0, endpoint),
                            initN2=initN2[j],
                            rateN2=rateN2[j],
                            ppN2=ppN2,
                            initHe=initHe[j],
                            rateHe=rateHe[j],
                            ppHe=ppHe,
                            M0=M0[j])$root
        }
      }
      result[i] <- min(val)
    }
    return(result)
  }
}

