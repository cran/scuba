#
# 	dive.R
#
#	$Revision: 1.18 $	$Date: 2008/06/20 21:54:29 $
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
dive <- function(..., begin=0, end=0, tanklist=NULL) {
  X <- list(...)
  stopifnot(is.numeric(begin) && length(begin) == 1)
  stopifnot(is.numeric(end) && length(end) == 1)
  tanks.given <- !is.null(tanklist)

  # initialise current state and rules
  state <- data.frame(time=0, depth=begin,
                      fO2=air$fO2, fN2=air$fN2, fHe=air$fHe, tankid=1)
  rate.down <- descent(30)  # default descent rate
  rate.up <- ascent(18)    # default ascent rate
  if(!tanks.given)
    tanklist <- list(air)

  # set up data frame to record the entire dive profile
  # NB: gas fractions apply to the time interval [i,i+1]
  df <- state
  now <- 1

  # Some entries of X may be labelled (tag=value)
  tagged <- nzchar(names(X))
  if(length(tagged) == 0) tagged <- rep(FALSE, length(X)) 
  
  change <- function(state, duration=NULL, depth=NULL,
                     fO2=NULL, fN2=NULL, fHe=NULL, tankid=NULL) {
    if(!is.null(duration)) state$time <- state$time + duration
    if(!is.null(depth))    state$depth <- depth
    if(!is.null(fO2))      state$fO2 <- fO2
    if(!is.null(fN2))      state$fN2 <- fN2
    if(!is.null(fHe))      state$fHe <- fHe
    if(!is.null(tankid))   state$tankid <- tankid
    return(state)
  }

  gasequal <- function(g1, g2) { identical(all.equal(g1, g2), TRUE) }

  # examine arguments
  for(i in seq(X)) {
    Y <- X[[i]]
    if(is.rate(Y)) # change the ascent or descent rate
      if(Y$up) rate.up <- Y else rate.down <- Y
    else if(is.gas(Y)) {
      # switch to new gas for next time interval
      if(tanks.given)
        stop(paste("If tanklist is given, gas switches should indicate",
                   "which tank is used (tank=number or tank=name)"))
      matched <- unlist(lapply(tanklist, gasequal, g2=Y))
      if(any(matched)) 
        tankid <- min(which(matched))
      else {
        # new gas
        if(now == 1) 
          tanklist <- list(Y)
        else 
          tanklist <- append(tanklist, list(Y))
        tankid <- length(tanklist)
      }
      state <- change(state, fO2=Y$fO2, fN2=Y$fN2, fHe=Y$fHe, tankid=tankid)
      df[now, ] <- state
    } else if(tanks.given && tagged[i] && names(X)[i] == "tank") {
      # switch to new tank Y for next time interval
      if(is.null(tanklist))
        stop("Cannot change tank: no tanks specified")
      Z <- tanklist[[Y]]
      state <- change(state, fO2=Z$fO2, fN2=Z$fN2, fHe=Z$fHe, tankid=Y)
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

      newdf <- data.frame(time=newtimes,
                          depth=newdepths,
                          fO2=state$fO2,
                          fN2=state$fN2,
                          fHe=state$fHe)
      if(!is.null(tanklist))
        newdf <- cbind(newdf, tankid=state$tankid)
      df <- rbind(df, newdf)
      now <- now + nsteps
      state <- df[now, ]
    } else if(is.dive(Y)) {
      #### dive object
      Ydata <- Y$data
      # First ascend/descend to depth (on current gas)
      newdepth <- Ydata$depth[1]
      if(newdepth != state$depth) {
        duration <- timetaken(state$depth, newdepth, rate.up, rate.down)
        state <- change(state, duration, newdepth)
        df <- rbind(df, state)
        now <- now + 1
      }
      # Now reconcile breathing gases.
      if(now == 1 && !tanks.given) {
        # we haven't started diving yet, so no tanks have been breathed.
        # Take tank list from Y and use gas specified by Y
        tanklist <- Y$tanklist
        tanks.given <- Y$tanks.given
      } else {
        # Concatenate tank lists and convert tank ID's to integers
        tid <- as.integer(df$tankid)
        Ytid <- as.integer(Ydata$tankid)
        tanklist <- append(tanklist, Y$tanklist)
        df$tankid <- tid
        Ydata$tankid <- Ytid + max(tid)
      }
      # Then tack on the new dive
      Ydata$time <- (Ydata$time - Ydata$time[1]) + state$time
      df <- rbind(df, Ydata)
      state <- df[nrow(df), ]
      now <- nrow(df) + 1
    } else 
    stop(paste("The format of argument", i, "is not recognised\n"))
  }
  if(state$depth != end) {
    # Don't forget to surface.. (or whatever)
    duration <- timetaken(state$depth, end, rate.up, rate.down)
    state <- change(state, duration, end)
    df <- rbind(df, state)
    now <- now + 1
  }
  # was it an air dive?
  airdive <- all(with(df, gasnames(fO2, fN2, fHe)) == "air")

  # dive depth/time profile
  pro <- approxfun(df$time, df$depth, method="linear", rule=2, ties="ordered")
  # 
  result <- list(profile=pro, data=df, airdive=airdive, tanklist=tanklist,
                 tanks.given=tanks.given)
  class(result) <- c("dive", class(result))
  return(result)
}

is.dive <- function(x) { inherits(x,"dive")}

depths.dive <- function(d) {
  stopifnot(is.dive(d))
  d$data$depth
}

times.dive <- function(d) {
  stopifnot(is.dive(d))
  d$data$time
}

durations.dive <- function(d) {
  stopifnot(is.dive(d))
  diff(times.dive(d))
}

"depths.dive<-" <- function(d, value) {
  stopifnot(is.dive(d))
  stopifnot(is.vector(value) && is.numeric(value))
  if(length(value) != nrow(d$data))
    stop("incorrect length of replacement vector for depths")
  d$data$depth <- value
  # update dive depth/time profile
  d$profilepro <- approxfun(d$data$time, d$data$depth,
                            method="linear", rule=2, ties="ordered")
  # 
  return(d)
}

"times.dive<-" <- function(d, value) {
  stopifnot(is.dive(d))
  stopifnot(is.vector(value) && is.numeric(value))
  if(length(value) != nrow(d$data))
    stop("incorrect length of replacement vector for times")
  d$data$time <- value
  # update dive depth/time profile
  d$profilepro <- approxfun(d$data$time, d$data$depth,
                            method="linear", rule=2, ties="ordered")
  return(d)
}

"durations.dive<-" <- function(d, value) {
  stopifnot(is.dive(d))
  stopifnot(is.vector(value) && is.numeric(value))
  if(length(value) != nrow(d$data) - 1)
    stop("incorrect length of replacement vector for durations")
  d$data$time <- d$data$time[1] + cumsum(c(0, value))
  # update dive depth/time profile
  d$profilepro <- approxfun(d$data$time, d$data$depth,
                            method="linear", rule=2, ties="ordered")
  return(d)
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
    labels <- with(x$data, gasnames(fO2, fN2, fHe))[-nrow(x$data)]
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
  if(any(iii <- is.infinite(times)))
    watchtimes[iii] <- "Inf"
  
  if(x$airdive) {
    # air dive: print time and depth profiles
    cat("gas: air\n")
    print(data.frame(time=watchtimes, depth=depths))
  } else if(x$tanks.given) {
    # dive using specified tanks
    print(data.frame(time=watchtimes, depth=depths, tank=factor(x$data$tankid),
                     stringsAsFactors=FALSE), quote=FALSE)
    cat("\nTank contents:\n")
    tk <- x$tanklist
    gasblurbs <- unlist(lapply(tk, as.character, full=TRUE))
    tanknames <- names(tk)
    if(!any(nzchar(tanknames))) tanknames <- paste(1:length(tk))
    for(i in 1:length(gasblurbs))
      cat(paste(tanknames[i], ":\t", gasblurbs[i], "\n", sep=""))
  } else {
    # gases specified along the way.
    gases <- with(x$data, gasnames(fO2, fN2, fHe))
    if(length(unique(gases)) == 1) {
      cat(paste("gas:", gases[1], "\n"))
      print(data.frame(time=watchtimes, depth=depths))
    } else
    print(data.frame(time=watchtimes,
                     depth=depths,
                     gas=factor(gases),
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
  gases <- with(object$data, gasnames(fO2, fN2, fHe))
  flat <- (diff(depths) == 0 & durations >= 1)
  stages <- data.frame(depth=(depths[-n])[flat],
                       time=durations[flat],
                       gas=(gases[-n])[flat])
  z <- list(depths=depths,
            times=times,
            totaltime=totaltime,
            stages=stages,
            maxdepth=max(depths),
            meandepth=meandepth,
            airdive=object$airdive,
            gasnames=gases)
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


#################################################################
#
# Manipulating dives
#

chop.dive <- function(d, t0=0, t1=max(times.dive(d))) {
  stopifnot(is.dive(d))
  stopifnot(t1 > t0)
  times <- times.dive(d)
  df <- d$data
  if(t0 < 0)
    stop("The start time t0 is before the beginning of the dive")
  if(t1 > max(times))
    stop("The termination time t1 is later than the end of the dive")

  intervening <- (times >= t0 & times <= t1)
  middle <- any(intervening)
  if(middle) {
    # extract data from waypoints between t0 and t1 inclusive
    ra <- range(which(intervening))
    first <- ra[1]
    last <- ra[2]
    dfOUT <- df[first:last, ]
  } else {
    # there are no intervening waypoints
    dfOUT <- NULL
    # index for data that apply to dive
    first <- last <- max(which(times <= t0))
  }
  
  if(times[first] != t0) {
    # start time is between two waypoints
    # interpolate to find depth
    dep0 <- d$profile(t0)
    # tack on an initial row
    df0 <- df[first,]
    df0$time <- t0
    df0$depth <- dep0
    dfOUT <- rbind(df0, dfOUT)
  }
  if(times[last] != t1) {
    # finish time is between two waypoints
    # interpolate to find depth
    dep1 <- d$profile(t1)
    # tack on a final row
    df1 <- df[last,]
    df1$time <- t1
    df1$depth <- dep1
    dfOUT <- rbind(dfOUT, df1)
  }
  # reset the clock
  dfOUT$time <- with(dfOUT, time - time[1])
  # create dive profile 
  pro <- approxfun(dfOUT$time, dfOUT$depth,
                   method="linear", rule=2, ties="ordered")
  airdive <- all(with(dfOUT, gasnames(fO2, fN2, fHe)) == "air")
  result <- d
  result$profile <- pro
  result$data    <- dfOUT
  result$airdive <- airdive
  return(result)
}

dive.segment <- function(d, i) {
  # make a `dive' out of one segment (between two waypoints) of a dive d
  df <- d$data[c(i, i+1), ]
  airdive <- all(with(df, gasnames(fO2, fN2, fHe)) == "air")
  pro <- approxfun(df$time, df$depth, method="linear", rule=2, ties="ordered")
  result <- d
  result$profile <- pro
  result$data    <- df
  result$airdive <- airdive
  return(result)
}


###########################################################
# tank lists

tanklist <- function(d) {
  stopifnot(is.dive(d))
  return(d$tanklist)
}

"tanklist<-" <- function(d, value) {
  if(!is.dive(d))
    stop("In tanklist(d) <- value, d must be a dive object")
  if(!is.list(value) || !all(unlist(lapply(value, is.gas))))
    stop("In tanklist(d) <- value, value must be a list of gases")
  tk <- d$tanklist
  if(is.null(tk))
    stop("Dive does not contain any tank information")
  if(length(tk) != length(value))
    stop("Lengths of new and old tanklists do not match")

  # reconcile names of tanks
  hasoldnames <- any(nzchar(names(tk)))
  hasnewnames <- any(nzchar(names(value)))
  if(!hasnewnames && hasoldnames)
    names(value) <- names(tk)

  d$tanks.given <- TRUE
  d$tanklist <- value
  df <- d$data
  for(i in 1:nrow(df)) {
    id <- df$tankid[i]
    g <- value[[id]]
    df$fO2[i] <- g$fO2
    df$fN2[i] <- g$fN2
    df$fHe[i] <- g$fHe
  }

  if(hasnewnames) 
    df$tankid <- factor(names(value)[as.integer(df$tankid)],
                        levels=names(value))

  d$data <- df
  d$airdive <- all(unlist(lapply(value, is.air)))
  return(d)
}

allspecies <- function(d, inert=TRUE) {
  stopifnot(is.dive(d))
  gf <- d$data[, c("fO2", "fN2", "fHe")]
  present <- (apply(as.matrix(gf), 2, max) > 0)
  spec <- c("O2", "N2", "He")
  out <- spec[as.logical(present)]
  if(inert)
    out <- out[out != "O2"]
  return(out)
}
