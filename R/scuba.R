#
# 	scuba.R
#
#	$Revision: 1.14 $	$Date: 2007/04/16 07:08:08 $
#
###################################################################
#  
# Dive profiles
#
# e.g.
#      dive(c(21,30),c(5,5),c(0,120),c(18,30),c(5,5))
#	  means: 21 m for 30 min, 5@5, 2 hour surface, 18 m for 30 min, 5@5
#
#      dive(c(25,30,0.32))  means EANx 32
#
dive <- function(...) {
  X <- list(...)
  # vectors describing the dive
  times <- 0     # start at time zero
  depths <- 0    # start at surface
  current <- 1
  gases <- list() # gas for first interval not yet determined
  # current values
  gas.current <- air
  rate.down <- descent(30)  # default descent rate
  rate.up <- ascent(18)    # default ascent rate
  # examine arguments
  for(i in seq(X)) {
    Y <- X[[i]]
    if(is.rate(Y)) # change the ascent or descent rate
      if(Y$up) rate.up <- Y else rate.down <- Y
    else if(is.gas(Y)) # switch to new gas
      gas.current <- Y
    else if(is.vector(Y)) {
      newdepth <- Y[1]
      # First ascend or descend to this depth
      if(newdepth != depths[current]) {
        duration <- timetaken(depths[current], newdepth, rate.up, rate.down)
        gases[[current]] <- gas.current
        times <- c(times, times[current] + duration)
        depths <- c(depths, newdepth)
        current <- current + 1
      }
      # Now stay at this depth for the specified time
      if(length(Y) > 1) {
        duration <- Y[2]
        gases[[current]] <- gas.current
        times <- c(times, times[current] + duration)
        depths <- c(depths, depths[current])
        current <- current + 1
      }
    }
  }
  if(depths[current] > 0) {
    # Don't forget to surface...
    gases[[current]] <- gas.current
    duration <- timetaken(depths[current], 0, rate.up, rate.down)
    times <- c(times, times[current] + duration)
    depths <- c(depths, 0)
    current <- current + 1
  }
  # was it an air dive?
  airdive <- all(unlist(lapply(gases, is.air)))
  
  # dive depth/time profile
  pro <- approxfun(times, depths, method="linear", rule=2, ties="ordered")
  # 
  result <- list(profile=pro, gases=gases, airdive=airdive)
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
  
plot.dive <- function(x, ..., main=deparse(substitute(x)), verticals=TRUE) {
  times <- times.dive(x)
  depths <- depths.dive(x)
  plot(times, -depths, type="o", main=main,
       axes=FALSE,
       xlab="Time (minutes)", ylab="Depth (metres)", lwd=2, ...)
  axis(1, at=pretty(range(times)))
  yp <- pretty(range(depths))
  axis(2, at=-yp, labels=paste(yp))

  if(!x$airdive) {
    isair <- unlist(lapply(x$gases, is.air))
    fO2 <- unlist(lapply(x$gases, function(g) { g$fO2 }))
    labels <- ifelse(isair, "air", paste("EAN ", round(100 * fO2)))
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
                      ifelse(secs== 0, "00", paste(secs)), sep="")
  } else 
    watchtimes <- paste(round(times))
  
  if(x$airdive) {
    cat("gas: air\n")
    print(data.frame(time=watchtimes, depth=depths))
  } else {
    gasnames <- unlist(lapply(x$gases, function(g) { summary(g)$name }))
    print(data.frame(time=watchtimes, depth=paste(depths), gas=c(gasnames, "")),
          quote=FALSE)
  }
  invisible(NULL)
}

summary.dive <- function(object, ...) {
  depths <- depths.dive(object)
  times <- times.dive(object)
  n <- length(depths)
  totaltime <- times[n]
  durations <- diff(times)
  middepths <- (depths[-1] + depths[-n])/2
  meandepth <- sum(durations * middepths)/totaltime
  gasnames <- unlist(lapply(object$gases, function(g) { summary(g)$name }))
  flat <- (diff(depths) == 0)
  stages <- data.frame(depth=depths[flat],
                       time=durations[flat],
                       gas=gasnames[flat])
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
  else 
    cat(paste("\nGases:", paste(unique(x$gasnames), collapse=", "), "\n"))
  cat(paste("Total dive time:", round(x$totaltime,1), "minutes\n"))
  cat("Stages:\n")
  if(x$airdive)
    print(x$stages[,1:2])
  else
    print(x$stages)
  cat(paste("\nMean depth", round(x$meandepth,1), "metres\n"))
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

print.gas <- function(x, ...) {
  name <- if(is.air(x)) "air" else if(x$fO2 == 1) "100% O2" else paste("EAN", 100 * x$fO2, "%")
  cat(paste(name, "\n"))
  invisible(NULL)
}

summary.gas <- function(object, ...) {
  fo <- object$fO2
  na <- if(is.air(object)) "air" else if(fo == 1) "100% O2" else paste("EAN", 100 * fo, "%")
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
	 prevstate=rep(0.79,length(halftimes))) 
{
  stopifnot(is.dive(d))
  if(length(prevstate) != length(halftimes))
    stop("the length of \'prevstate\' does not match the number of compartments in the model")
  ratecoefs <- log(2)/halftimes

  # 'state' is the vector of nitrogen tensions (ata) in each compartment
  state <- prevstate

  times <- times.dive(d)
  depths <- depths.dive(d)
  fN2 <- unlist(lapply(d$gases, function(g) { 1 - g$fO2 }))
  durations <- diff(times)
  n <- length(times)

  for(i in 1:(n-1)) {
    d0 <- depths[i]
    d1 <- depths[i+1]
    tim <- durations[i]
    k1 <- fN2[i] * (d0/10 + 1)
    k2 <- fN2[i] * ((d1-d0)/10) / tim
    f <- exp(- ratecoefs * tim)
    state <- state * f + (k1 - k2/ratecoefs) * (1 - f) + k2 * tim
  }
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
function(d)
{
  times <- times.dive(d)
  depths <- depths.dive(d)
  n <- length(depths)

  fO2 <- unlist(lapply(d$gases, function(g) { g$fO2 }))
  fO2 <- c(fO2, fO2[n-1])
  pO2 <- fO2 * (depths/10 + 1)
  if(any(pO2 > 1.6))
    warning("O2 partial pressure exceeded 1.6")
  else if(any(pO2 > 1.4))
    warning("O2 partial pressure exceeded 1.4")
  toxic <- (pO2 > 0.5)
  if(!any(toxic)) return(0)

  # calculations for each interval
  durations <- diff(times)
  flat <- (diff(depths) == 0)
  tox <- toxic[-n] | toxic[-1]
  powerit <- function(x) { exp(0.83 * log(x)) }

  otu.sum <- 0
  for(i in seq(tox)[tox]) {
    dura <- durations[i]
    p0 <- pO2[i]
    p1 <- pO2[i+1]
    if(flat[i])
      # integrate the constant (2 *(p02-0.5))^0.83
      otu <- dura * powerit(2 * (p0 - 0.5))
    else if(toxic[i] && toxic[i+1]) {
      # integrate (a + bx)^0.83)
      a <- 2 * p0 - 1
      b <- 2 * (p1-p0)/dura
      otu <- (1/(1.83 * b)) * (powerit(a + b * dura) - powerit(a))
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
      otu <- (1/(1.83 * b)) * (powerit(a + b * dura) - powerit(a))
    }
    otu.sum <- otu.sum + otu
  }
  return(otu.sum)
}


#################################################################
#
# Manipulating dives
#

chop.dive <- function(d, tim) {
  stopifnot(is.dive(d))
  stopifnot(tim > 0)
  times <- times.dive(d)
  depths <- depths.dive(d)
  gases <- d$gases
  if(tim > max(times))
    stop("The interruption time \"tim\" is later than the end of the dive")
  last <- max(seq(times)[times <= tim])
  if(times[last] == tim) {
    if(last==1)
      stop("Zero-duration dive")
    else {
      times <- times[1:last]
      depths <- depths[1:last]
      gases <- gases[1:(last-1)]
    }
  } else {
      times <- c(times[1:last], tim)
      depths <- depths[c(1:last, last)]
      gases <- gases[c(1:(last-1), last-1)]
  }
  pro <- approxfun(times, depths, method="linear", rule=2, ties="ordered")
  airdive <- all(unlist(lapply(gases, is.air)))
  result <- list(profile=pro, gases=gases, airdive=airdive)
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

  # determine (approx) maximum relative saturation over entire dive
  # so that barplot scale is fixed
  
  maxrelsat <- 1
  for(ti in times.dive(d))
    if(ti > 0) {
      div <- chop.dive(d, ti)
      hal <- haldane(div, halftimes=HalfT)
      rel <- hal/Mvals
      maxrelsat <- max(max(rel), maxrelsat)
    }

  # loop ..
  
  while(length(xy <- locator(1)) > 0) {
    tim <- xy$x
    if(tim <= 0) 
      cat("time < 0; pick again\n")
    else {
      tim <- min(tim, max(times.dive(d)))
      div <- chop.dive(d, tim)
      hal <- haldane(div, halftimes=HalfT)
      rel <- hal/Mvals
      otu <- oxtox(div)
      barplot(rel,
              xlab=paste("Tissues (", model, ")", sep=""),
              ylab="Relative saturation",
              main=paste("Time=", round(tim), "min\n",
                          "Accumulated oxygen toxicity", round(otu, 1)),
              ylim=range(c(0, maxrelsat)),
              names.arg=c("Fast", rep("", length(rel)-2), "Slow"))
      abline(h=1, lty=3, col="red")
      plot(d)
      abline(v=tim, lty=2, col="green")
    }
  }
  par(oldpar)
  hal
}

