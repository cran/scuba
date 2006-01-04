\name{dive}
\alias{dive}
\title{Define a Dive Profile}
\description{
  Define a dive profile.
}
\usage{
  dive(\dots)
}
\arguments{
  \item{\dots}{
    Any number of arguments, specifying a sequence of events
    that make up the dive. The arguments may specify depths,
    time spent at each depth, ascent or descent rates,
    and gas switches. See Details.
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
      created by the function \code{\link{nitrox}}.
      The diver switches to this gas.
    }
    \item{an ascent/descent rate object}{
      created by the functions \code{\link{ascent}} or
      \code{\link{descent}}. This resets the default
      rate of ascent or descent.
    }
  }
  Dives are assumed to start and finish at the surface.

  Initially the descent rate is set to 30 metres per minute,
  the descent rate is 18 metres per minute, and the breathing gas is
  air. These settings may be changed during the dive
  by the \code{...} arguments.

  Once an object of class \code{"dive"} has been created,
  it can be plotted and printed (using \code{\link{plot.dive}}
  and \code{\link{print.dive}}). The nitrogen saturation can be
  computed using \code{\link{haldane}} and the
  cumulative oxygen toxicity using \code{\link{oxtox}}.
}
\section{Warnings}{
  Not suitable for representing altitude dives. Not applicable to
  gas mixtures other than air and nitrox (oxygen-nitrogen mixtures).
}
\seealso{
  \code{\link{nitrox}},
  \code{\link{ascent}},
  \code{\link{descent}},
  \code{\link{haldane}},
  \code{\link{plot.dive}}
}
\examples{
  # Dive to 25 m for 20 min followed by safety stop at 5 metres for 3 min
  d <- dive(c(25,20),c(5,3))
  plot(d)

  # Bounce dive to 20 metres
  d <- dive(20)

  # GAS SWITCHING
  # Dive to 18 m for 30 min on air,
  # switch to Nitrox 36, ascend to 5 metres, safety stop
  d <- dive(c(18, 30), nitrox(0.36), c(5,3))
  # Same as above, but ascend to 5 m on air, then switch gas
  d <- dive(c(18, 30), 5, nitrox(0.36), c(5,3))

  # ASCENT RATES
  # Ascent rate 18 m/min below 9 metres, 6m/min above 9 metres
  d <- dive(c(30, 12), ascent(18), 9, ascent(6), c(5,3))
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
