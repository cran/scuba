\name{dive}
\alias{dive}
\title{Define a Dive Profile}
\description{
  Define a dive profile.
}
\usage{
  dive(\dots, begin=0, end=0, tanklist=NULL)
}
\arguments{
  \item{\dots}{
    Any number of arguments, specifying a sequence of events
    that make up the dive. The arguments may specify depths,
    time spent at each depth, ascent or descent rates,
    and gas switches. See Details.
  }
  \item{begin, end}{
    The depths at the start and finish of the dive (normally zero
    indicating that the dive starts and finishes at the surface).
  }
  \item{tanklist}{
    Optional list of the gases contained in each tank,
    for a dive with multiple tanks.
  }
}
\value{
  An object of class \code{"dive"} describing the dive profile
  and the breathing gases used.
}
\details{
  This function creates 
  an object of class \code{"dive"} which represents a scuba dive.
  The dive profile is assumed to be piecewise linear, that is,
  the graph of depth against time is a broken straight line. 
  Dives are assumed to start and finish at the surface.

  The arguments \code{...} specify a succession of events
  that make up the dive. Each argument can be
  \describe{
    \item{a vector of length 2}{
      interpreted as \code{c(depth, duration)}
      where \code{depth} gives the depth in metres
      and \code{duration} the length of stay at this depth, in minutes.
    }
    \item{a single number}{
      interpreted as a waypoint depth in metres. 
      The diver will ascend or descend to this depth, at the
      current default rate of ascent or descent.
    }
    \item{an object of class \code{"gas"}}{
      specifying a breathing gas. Such objects are
      created by the functions \code{\link{nitrox}} and \code{\link{trimix}}.
      The diver switches to this gas.
    }
    \item{an argument of the form \code{tank=n}}{specifying
      a switch to another tank or cylinder of breathing gas.
      This is available only when the tanks used in the dive
      have been specified by the argument \code{tanklist}.
      The value \code{n} should be a valid index for this list
      (either a serial number in the tank list or a character string
      matching one of the names in the tank list).
    }
    \item{an ascent/descent rate object}{
      created by the functions \code{\link{ascent}} or
      \code{\link{descent}}. This resets the default
      rate of ascent or descent.
    }
    \item{a data frame with 2 columns}{
      containing a dive profile, usually uploaded from a dive computer. The
      first column should contain the \emph{elapsed} times, and the second
      column contains the depths (in metres) measured at these times.
      The column of elapsed times can be either a character vector
      containing times in minutes-and-seconds format \code{mm:ss},
      or an integer vector containing elapsed times in \bold{seconds},
      or a vector of class \code{\link{difftime}}
      representing elapsed times in any time unit.
    }
    \item{an object of class \code{"dive"}}{
      representing a dive profile. This allows the user to paste
      dive profiles together.
    }
  }
  Dives are assumed to start and finish at the surface (this can be changed
  by specifying values for \code{begin} and \code{end}).

  Initially the descent rate is set to 30 metres per minute,
  the ascent rate is 18 metres per minute, and the breathing gas is
  air. These settings may be changed during the dive
  by the \code{...} arguments.

  A dive object may include
  periods spent at the surface (depth zero) and may therefore represent
  a succession of dives separated by surface intervals. 

  Once an object of class \code{"dive"} has been created,
  it can be plotted and printed (using \code{\link{plot.dive}}
  and \code{\link{print.dive}}). The nitrogen saturation can be
  computed using \code{\link{haldane}} and the
  cumulative oxygen toxicity using \code{\link{oxtox}}.
}
\section{Modifying a dive object}{
  The depths and elapsed times at each waypoint during the dive
  can be extracted and changed using \code{\link{depths.dive}}
  and \code{\link{times.dive}}. 

  It is possible to alter the \emph{breathing gases} used in a dive \code{d},
  yielding a new dive object. 
  You simply modify the tanklist using \code{\link{tanklist}(d) <- value}. 
  This makes it possible to study the effect of conducting the same
  dive profile with a different breathing gas.
}
\section{Warnings}{
  Not suitable for representing altitude dives. 
}
\seealso{
  \code{\link{nitrox}},
  \code{\link{ascent}},
  \code{\link{descent}},
  \code{\link{haldane}},
  \code{\link{plot.dive}},
  \code{\link{times.dive}},
  \code{\link{durations.dive}},
  \code{\link{depths.dive}},
  \code{\link{chop.dive}},
  \code{\link{tanklist}}.
}
\examples{
  # Dive to 25 m for 20 min followed by safety stop at 5 metres for 3 min
  d <- dive(c(25,20),c(5,3))
  plot(d)

  # Bounce dive to 20 metres
  d <- dive(20)

  # Two dives separated by a one-hour surface interval
  d <- dive(c(30,15),c(9,2),c(5,5),c(0,60),c(12,60),c(5,5))

  # ASCENT RATES
  # Ascent rate 18 m/min below 9 metres, 6m/min above 9 metres
  d <- dive(c(30, 12), ascent(18), 9, ascent(6), c(5,3))

  # UPLOADED DIVE PROFILE
  data(baron)
  pro <- baron[, 1:2]
  d <- dive(pro)
  plot(d)

  # GAS USE
  # 30-metre dive on Nitrox 32
  d <- dive(nitrox(0.32), c(30,20), c(5,5))

  # GAS SWITCHING
  # Dive to 18 m for 30 min on air,
  # switch to Nitrox 36, ascend to 5 metres, safety stop
  d <- dive(c(18, 30), nitrox(0.36), c(5,3))
  # Same as above, but ascend to 5 m on air, then switch gas
  d <- dive(c(18, 30), 5, nitrox(0.36), c(5,3))

  # SWITCHING TO SPECIFIC TANKS
  d <- dive(tanklist=list(main=air, deco=nitrox(0.50)),
            tank="main", c(30, 20), 5, tank="deco", c(5,10))

  # Descend to 5 metres on pure oxygen, switch to Trimix,
  # descend to 30 metres, remain 40 minutes, ascend to 6 metres,
  # switch to pure oxygen, ascend to 5 metres, decompress 10 minutes,
  # surface and continue breathing pure oxygen for 10 minutes

  d <- dive(tanklist=list(travel=trimix(0.15, 0.5), deco=nitrox(1)),
            tank="deco", 5, tank="travel", c(30,40), 6, tank="deco",
            c(5,10), c(0,10))

}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
