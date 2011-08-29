\name{durations.dive}
\alias{durations.dive}
\alias{durations.dive<-}
\title{Durations of time between each waypoint of a dive}
\description{
  Extracts, or alters, the durations of each interval
  between waypoints during a dive.
}
\usage{
durations.dive(d)
durations.dive(d) <- value
}
\arguments{
  \item{d}{A dive (object of class \code{"dive"}).}
  \item{value}{A numeric vector containing durations in minutes.}
}
\details{
  An object of class \code{"dive"} represents a scuba dive. It is
  created by the function \code{\link{dive}}. A dive is defined as a
  series of waypoints occurring at specified depths and times.

  The duration of the time interval between each successive pair of
  waypoints is returned by \code{durations.dive}. The
  assignment \code{durations.dive(d) <- value} alters these durations,
  provided the new vector \code{value} has the same length as the old one.
}
\value{
  \code{durations.dive} returns a numeric vector containing the
  durations of intervals between each waypoint, in minutes.
}
\seealso{
  \code{\link{depths.dive}},
  \code{\link{times.dive}},
  \code{\link{dive}}.
}
\examples{
   d <- dive(c(30,20), c(5,5))
   d
   durations.dive(d)
   # what if we stayed at the bottom for 25 minutes instead of 20?
   durations.dive(d)[2] <- 25
   d
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
