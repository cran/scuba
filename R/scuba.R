#
# 	scuba.R
#
#	$Revision: 1.26 $	$Date: 2007/08/29 18:30:04 $
#
###################################################################
#  
# Dive profiles
#
# e.g.
#      dive(c(21,30),c(5,5),c(0,120),c(18,30),c(5,5))
#	  means: 21 m for 30 min, 5@5, 2 hour surface, 18 m for 30 min, 5@5
#
#      dive(nitrox(0.32), c(25,30), 5, nitrox(0.5), c(5,5))
#                 means dive on EANx 32 to 25 metres for 30 min,
#                 ascending to 5 metres, switching to EANx 50,
#                 5@5 then surfacing.
#
dive <- function(...) {
  X <- list(...)
  
  # initialise current state and rules
  state <- data.frame(time=0, depth=0, fO2=air$fO2, fN2=air$fN2)
  rate.down <- descent(30)  # default descent rate
  rate.up <- ascent(18)    # default ascent rate
  
  # set up data frame to record the entire dive profile
  # NB: gas fractions apply to the time interval [i,i+1]
  df <- state
  now <- 1

  change <- function(state, duration=NULL, depth=NULL, fO2=NULL, fN2=NULL) {
    if(!is.null(duration)) state$time <- state$time + duration
    if(!is.null(depth))    state$depth <- depth
    if(!is.null(fO2))      state$fO2 <- fO2
    if(!is.null(fN2))      state$fN2 <- fN2
    return(state)
  }
  
  # examine arguments
  for(i in seq(X)) {
    Y <- X[[i]]
    if(is.rate(Y)) # change the ascent or descent rate
      if(Y$up) rate.up <- Y else rate.down <- Y
    else if(is.gas(Y)) { # switch to new gas for next time interval
      state <- change(state, fO2=Y$fO2, fN2=Y$fN2)
      df[now, ] <- state
    } else if(is.vector(Y)) {
      if(length(Y) > 2)
        stop("Vector of length > 2 is not a recognised format")
      newdepth <- Y[1]
      # First ascend or descend to this depth (making new waypoint)
      if(newdepth != state$depth) {
        duration <- timetaken(state$depth, newdepth, rate.up, rate.down)
        state <- change(state, duration, newdepth)
        df <- rbind(df, state)
        now <- now + 1
      }
      # Now stay at this depth for the specified time (making new waypoint)
      if(length(Y) > 1) {
        duration <- Y[2]
        state <- change(state, duration)
        df <- rbind(df, state)
        now <- now + 1
      }
    } else if(is.data.frame(Y)) {
      if(ncol(Y) != 2)
        stop("Data frame should have 2 columns")
      newtimes <- Y[,1]
      newdepths <- Y[,2]
      nsteps <- nrow(Y)

      if(inherits(newtimes, "difftime")) {
        # convert to minutes
        units(newtimes) <- "mins"
        newtimes <- as.numeric(newtimes)
      } else if(is.numeric(newtimes)) {
        # assume seconds -- convert to minutes
        message("Elapsed times are assumed to be given in seconds")
        newtimes <- newtimes/60
      } else if(is.character(newtimes)) {
        # assume mm:ss format
        newtimes <- as.difftime(newtimes, format="%M:%S")
        units(newtimes) <- "mins"
        newtimes <- as.numeric(newtimes)
      } else 
        stop("Unknown format for elapsed times")
      
      if(!all(diff(newtimes) > 0))
        stop("Elapsed times are not an increasing sequence")

      newdf <- data.frame(depth=newdepths,
                          time=newtimes,
                          fO2=state$fO2,
                          fN2=state$fN2)
      df <- rbind(df, newdf)
      now <- now + nsteps
      state <- df[now, ]
    } else
    stop(paste("The format of argument", i, "is not recognised\n"))
  }
  if(state$depth > 0) {
    # Don't forget to surface...
    duration <- timetaken(state$depth, 0, rate.up, rate.down)
    state <- change(state, duration, 0)
    df <- rbind(df, state)
    now <- now + 1
  }
  # was it an air dive?
  airdive <- all(df$fO2 == air$fO2)

  # dive depth/time profile
  pro <- approxfun(df$time, df$depth, method="linear", rule=2, ties="ordered")
  # 
  result <- list(profile=pro, data=df, airdive=airdive)
  class(result) <- c("dive", class(result))
  return(result)
}

is.dive <- function(x) { inherits(x,"dive")}

depths.dive <- function(d) {
  eval(expression(y), env=environment(d$profile))
}

times.dive <- function(d) {
  eval(expression(x), env=environment(d$profile))
}
  
plot.dive <- function(x, ...,
                      main=deparse(substitute(x)),
                      show.gases=TRUE, verticals=TRUE) {
  times <- times.dive(x)
  depths <- depths.dive(x)
  if(mean(diff(times)) <= 1) {
    show.gases <- FALSE
    plot.type <- "l"
  } else 
    plot.type <- "o"
  myplot <- function(x, y, main, ..., type=plot.type, axes=FALSE,
                     xlab="Time (minutes)", ylab="Depth (metres)", lwd=2) {
    plot(x, y, main=main, type=type, axes=axes,
         xlab=xlab, ylab=ylab, lwd=lwd, ...)
  }
  myplot(times, -depths, main=main, ...)
  axis(1, at=pretty(range(times)))
  yp <- pretty(range(depths))
  axis(2, at=-yp, labels=paste(yp))

  if(show.gases && !x$airdive) {
    fO2 <- x$data$fO2
    labels <- nitroxname(fO2)
    rates <- abs(diff(depths)/diff(times))
    aspect <- diff(range(depths))/diff(range(times))
    vertical <- (rates > 2 * aspect)
    n <- length(times)
    midx <- (times[-n] + times[-1])/2
    midy <- -x$profile(midx)
    text(midx[!vertical],  midy[!vertical], labels[!vertical], pos=3)
    if(verticals)
      text(midx[vertical], midy[vertical], labels[vertical], srt=90)
  }
  invisible(NULL)
}

print.dive <- function(x, ..., seconds=TRUE) {
  cat("Dive profile\n")
  depths <- depths.dive(x)
  times <- times.dive(x)
  if(seconds) {
    mins <- floor(times)
    secs <- floor((times * 60) %% 60)
    watchtimes <- paste(mins, ":",
                        ifelse(secs < 10, "0", ""),
                        secs, sep="")
  } else 
    watchtimes <- paste(round(times))
  
  if(x$airdive) {
    cat("gas: air\n")
    print(data.frame(time=watchtimes, depth=depths))
  } else {
    fO2 <- x$data$fO2
    if(length(unique(fO2)) == 1) {
      cat(paste("gas:", nitroxname(fO2[1]), "\n"))
      print(data.frame(time=watchtimes, depth=depths))
    } else
    print(data.frame(time=watchtimes,
                     depth=depths,
                     gas=nitroxname(fO2),
                     stringsAsFactors=FALSE), quote=FALSE)
  }
  return(invisible(NULL))
}

summary.dive <- function(object, ...) {
  depths <- depths.dive(object)
  times <- times.dive(object)
  n <- length(depths)
  totaltime <- times[n]
  durations <- diff(times)
  middepths <- (depths[-1] + depths[-n])/2
  meandepth <- sum(durations * middepths)/totaltime
  fO2 <- object$data$fO2
  gasnames <- nitroxname(fO2)
  flat <- (diff(depths) == 0 & durations >= 1)
  stages <- data.frame(depth=(depths[-n])[flat],
                       time=durations[flat],
                       gas=(gasnames[-n])[flat])
  z <- list(depths=depths,
            times=times,
            totaltime=totaltime,
            stages=stages,
            maxdepth=max(depths),
            meandepth=meandepth,
            airdive=object$airdive,
            gasnames=gasnames)
  class(z) <- c("summary.dive", class(z))
  return(z)
}
  
print.summary.dive <- function(x, ...) {
  cat(paste("Dive to", x$maxdepth, "metres"))
  if(x$airdive)
    cat(" on air\n")
  else {
    gases <- unique(x$gasnames)
    ngas <- length(gases)
    if(ngas <= 5) 
      cat(paste("\n", ngettext(ngas, "Gas: ", "Gases: "),
                paste(gases, collapse=", "), "\n", sep=""))
    else
      cat("More than 5 gas settings")
  }
  cat(paste("Total dive time:", round(x$totaltime,1), "minutes\n"))
  if(nrow(x$stages) > 0) {
    cat("Stages:\n")
    if(x$airdive)
      print(x$stages[,1:2])
    else
      print(x$stages)
  }
  cat(paste("Mean depth", round(x$meandepth,1), "metres\n"))
  return(invisible(NULL))
}
  
###################################################################
#  Gases
#

is.gas <- function(x) { inherits(x, "gas") }

nitrox <- function(fO2) {
  if(fO2 <= 0 || fO2 > 1)
    stop("fO2 should be a fraction between 0 and 1")
  g <- list(fO2=fO2, fN2=1-fO2)
  class(g) <- c("gas", class(g))
  return(g)
}

air <- nitrox(0.21)

is.air <- function(g) { g$fO2 == 0.21 }

nitroxname <- function(fO2) {
  ifelse(fO2 == 0.21, "air", ifelse(fO2 == 1, "100% O2", paste("EAN ", 100 * fO2, "%", sep="")))
}

print.gas <- function(x, ...) {
  name <- nitroxname(x$fO2)
  cat(paste(name, "\n"))
  invisible(NULL)
}

summary.gas <- function(object, ...) {
  fo <- object$fO2
  na <- nitroxname(fo)
  mo <- mod(object, ppO2max=1.4)
  mo <- round(mo, 1)
  z <- list(name=na, fO2=fo, fN2=1-fo, mod=mo)
  class(z) <- c("summary.gas", class(z))
  return(z)
}

print.summary.gas <- function(x, ...) {
  cat(paste(x$name, "\t(",
            100 * x$fO2, "% oxygen, ",
            100 * x$fN2, "% nitrogen)\n", sep=""))
  cat(paste("Maximum operating depth", x$mod, "metres\n"))
  return(invisible(NULL))
}

###################################################################
#  Ascent/descent rates/times
#

is.rate <- function(x) { inherits(x, "rate") }

rate <- function(speed=NULL, time=NULL, up=TRUE) {
  ngiven <- (!is.null(speed)) + !is.null(time)
  if(ngiven != 1)
    stop("Exactly one of \"speed\" and \"time\" should be specified")
  if(!is.null(speed) && speed <= 0)
    stop("speed should be a positive number (metres/minute)")
  if(!is.null(time) && time <= 0)
    stop("time should be a positive number (minutes)")
  ra <- list(speed=speed, time=time, up=up)
  class(ra) <- c("rate", class(ra))
  return(ra)
}

ascent <- function(speed=NULL, time=NULL) { return(rate(speed, time, TRUE)) }

descent <- function(speed=NULL, time=NULL) { return(rate(speed, time, FALSE)) }

print.rate <- function(x, ...) {
  cat(paste(
            if(x$up) "ascent" else "descent",
            if(!is.null(x$speed)) 
              paste("rate", x$speed, "metres/minute")
            else
              paste("time", x$time, "minutes"),
            "\n"))
  return(invisible(NULL))
}

timetaken <- function(start, finish, uprate, downrate) {
  if(!is.rate(uprate) || !uprate$up)
    stop("\"uprate\" should be an ascent rate")
  if(!is.rate(downrate) || downrate$up)
    stop("\"downrate\" should be a descent rate")
  rat <- if(start > finish) uprate else downrate
  if(!is.null(rat$time))
    return(rat$time)
  else
    return(abs(finish-start)/rat$speed)
}

        
####################################################################
#
#  Dive tables & models
#
"haldane" <- 
function(d,
	 halftimes={ data(Halftimes); Halftimes$DSAT}, 
	 prevstate=rep(0.79,length(halftimes)),
         progressive=FALSE) 
{
  stopifnot(is.dive(d))
  if(length(prevstate) != length(halftimes))
    stop("the length of \'prevstate\' does not match the number of compartments in the model")
  ratecoefs <- log(2)/halftimes

  # 'state' is the vector of nitrogen tensions (ata) in each compartment
  state <- prevstate

  times <- times.dive(d)
  depths <- depths.dive(d)
  fN2 <- d$data$fN2
  durations <- diff(times)
  n <- length(times)

  if(progressive) {
    profile <- as.data.frame(matrix(, nrow=n, ncol=length(halftimes)))
    profile[1,] <- prevstate
  }

  for(i in 1:(n-1)) {
    d0 <- depths[i]
    d1 <- depths[i+1]
    tim <- durations[i]
    if(tim > 0) {
      k1 <- fN2[i] * (d0/10 + 1)
      k2 <- fN2[i] * ((d1-d0)/10) / tim
      f <- exp(- ratecoefs * tim)
      state <- state * f + (k1 - k2/ratecoefs) * (1 - f) + k2 * tim
    }
    if(progressive)
      profile[i+1,] <- state
  }

  if(progressive)
    return(profile)
  else 
    return(state)
}

###################################################################
#  
# Nitrox 
#

"ead" <-
function(depth, g) {
  if(!is.gas(g))
    g <- nitrox(g)
  fO2 <- g$fO2
  (10 + depth) * (1- fO2)/0.79 - 10
}

"maxmix" <-
function(depth, ppO2max = 1.4) {
  fo <- pmin(1, ppO2max/(1 + depth/10))
  fo <- signif(fo, 2)
  nitrox(fo)
}

"mod" <-
function(g, ppO2max=1.4) {
  if(!is.gas(g))
    g <- nitrox(g)
  fO2 <- g$fO2
  10 * (ppO2max/fO2 - 1)
}

"eadtable" <-
function(g, ppO2max=1.4) {
  if(!is.gas(g))
    g <- nitrox(g)
  fO2 <- g$fO2
  depth <- 0:40
  EAD <- round(ead(depth, fO2), 1)
  EAD <- ifelse(EAD >= 0, as.character(EAD), " ")
  EAD <- ifelse(depth > mod(fO2, ppO2max), "Warning", EAD)
  data.frame(depth=depth,EAD=EAD)
}

"oxtox" <- 
function(d, progressive=FALSE)
{
  times <- times.dive(d)
  depths <- depths.dive(d)
  n <- length(depths)
  fO2 <- d$data$fO2

  # fO2[i] applies to interval (times[i], times[i+1])
  pO2start <- fO2 * (depths/10 + 1)
  pO2end   <- c(fO2[-n] * (depths[-1]/10 + 1), NA)
  pO2high  <- pmax(pO2start, pO2end, na.rm=TRUE)
  maxpO2 <- max(pO2high)
  if(maxpO2 > 1.6)
    warning("O2 partial pressure exceeded 1.6")
  else if(maxpO2 > 1.5)
    warning("O2 partial pressure exceeded 1.5")
  else if(maxpO2 > 1.4)
    warning("O2 partial pressure exceeded 1.4")

  # determine which intervals contain toxicity contributions
  toxicstart <- (pO2start > 0.5)
  toxicend   <- (pO2end   > 0.5)
  toxic <- (pO2high > 0.5)
  if(!any(toxic)) {
    if(!progressive) return(0) else return(rep(0, n))
  }

  # calculations for each interval
  durations <- diff(times)
  flat <- (diff(depths) == 0)
  powerit <- function(x) { exp(0.83 * log(x)) }
  Powerit <- function(x) { exp(1.83 * log(x)) }

  # vector of toxicity contributions
  dotu <- rep(0, n)
  
  # compute toxicity contributions for each toxic interval
  for(i in seq(n-1)[toxic[-n]]) {
    dura <- durations[i]
    p0 <- pO2start[i]
    p1 <- pO2end[i]
    if(flat[i])
      # integrate the constant (2 *(p02-0.5))^0.83
      dotu[i+1] <- dura * powerit(2 * (p0 - 0.5))
    else if(toxicstart[i] && toxicend[i]) {
      # entire interval is in the toxic zone
      # integrate (a + bx)^0.83)
      a <- 2 * p0 - 1
      b <- 2 * (p1-p0)/dura
      dotu[i+1] <- (1/(1.83 * b)) * (Powerit(a + b * dura) - Powerit(a))
    } else {
      # interval is only partly toxic
      # compute toxic subinterval
      t0 <- times[i]
      t1 <- times[i+1]
      tx <- t0 + (t1-t0) * (0.5-p0)/(p1-p0)
      # restrict to subinterval
      if(p0 < p1) {
        # [tx, t1]
        t0 <- tx
        p0 <- 0.5
        dura <- t1-tx
      } else {
        # [t0, tx]
        t1 <- tx
        p1 <- 0.5
        dura <- tx-t0
      }
      # integrate (a + bx)^0.83)
      a <- 2 * p0 - 1
      b <- 2 * (p1-p0)/dura
      dotu[i+1] <- (1/(1.83 * b)) * (Powerit(a + b * dura) - Powerit(a))
    }
  }
  if(progressive)
    return(cumsum(dotu))
  else
    return(sum(dotu))
}


#################################################################
#
# Manipulating dives
#

chop.dive <- function(d, tim) {
  stopifnot(is.dive(d))
  stopifnot(tim > 0)
  times <- times.dive(d)
  df <- d$data
  if(tim > max(times))
    stop("The interruption time is later than the end of the dive")
  last <- max(seq(times)[times <= tim])
  if(times[last] == tim) {
    # truncation time is an existing waypoint
    if(last==1)
      stop("Zero-duration dive")
    else 
      df <- df[1:last,]
  } else {
    # truncation time is between waypoints
    # interpolate to find depth
    dep <- d$profile(tim)
    df <- df[1:last,]
    df <- rbind(df, c(tim, dep, df$fO2[last], df$fN2[last]))
  }
  pro <- approxfun(df$time, df$depth, method="linear", rule=2, ties="ordered")
  airdive <- all(df$fO2 == 0.21)
  result <- list(profile=pro, data=df, airdive=airdive)
  class(result) <- c("dive", class(result))
  return(result)
}

###################################################################
#
#  Interactive display
#

showstates <- function(d, model="DSAT") {
  stopifnot(is.dive(d))
  stopifnot(is.character(model))

  model.table <- c("DSAT", "USN", "Workman", "ZH-L16A")
  k <- pmatch(model, model.table)
  if(is.na(k))
    stop(paste("Unrecognised model", sQuote(model)))
  model <- model.table[k]
  
  switch(model,
         DSAT={  
           data(Halftimes)
           data(Mvalues)
           HalfT <- Halftimes$DSAT
           Mvals <- Mvalues$ata$DSAT
         },
         USN={  
           data(Halftimes)
           data(Mvalues)
           HalfT <- Halftimes$USN
           Mvals <- Mvalues$ata$USN
         },
         Workman={
           data(Workman65)
           HalfT <- Workman65$halftime
           Mvals <- Workman65$M0/scuba.constants$msw.per.atm
         },
         "ZH-L16A"={
           data(BuehlmannL16A)
           HalfT <- BuehlmannL16A$halftime
           Mvals <- eval(expression(a + 1/b), envir=BuehlmannL16A)
         },
         stop(paste("Unrecognised model", sQuote(model)))
         )
         
  
  oldpar <- par(mfrow=c(1,2))
  frame()
  plot(d)
  timepoints <- times.dive(d)
  maxtime <- max(timepoints)
  ntissues <- length(HalfT)
  
  # precompute time profile of tissue saturation
  halprofile <- haldane(d, halftimes=HalfT, progressive=TRUE)

  # convert to relative saturation
  relhalprofile <- sweep(halprofile, 2, Mvals, "/")
  
  # determine (approx) maximum relative saturation over entire dive
  maxrelsat <- max(max(relhalprofile), 1)

  # precompute oxygen toxicity profile
  oxtoxprofile <- oxtox(d, progressive=TRUE)

  # event loop ..
  hal <- rep(NA, ntissues)

  while(length(xy <- locator(1)) > 0) {
    tim <- xy$x
    tim <- max(0, tim)
    tim <- min(tim, maxtime)
    ind <- min(which(timepoints >= tim))
    hal <- as.numeric(halprofile[ind,])
    rel <- as.numeric(relhalprofile[ind, ])
    otu <- oxtoxprofile[ind]
    
    pozzie <- barplot(rel,
                      xlab=paste("Tissues (", model, ")", sep=""),
                      ylab="Relative saturation",
                      main=paste("Time=", round(tim), "min\n",
                        "Accumulated oxygen toxicity", round(otu, 1)),
                      ylim=range(c(0, maxrelsat)),
                      names.arg=NULL)
    mtext(side=1, at=pozzie[1], line=1, text="Fast")
    mtext(side=1, at=pozzie[ntissues], line=1, text="Slow")
    abline(h=1, lty=3, col="red")
    plot(d)
    abline(v=tim, lty=2, col="green")
  }
  par(oldpar)
  hal
}

