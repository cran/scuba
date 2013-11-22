\name{times.dive}
\alias{times.dive}
\alias{times.dive<-}
\title{Elapsed times at each waypoint of a dive}
\description{
  Extracts, or alters, the elapsed time at each waypoint during a dive.
}
\usage{
times.dive(d)
times.dive(d) <- value
}
\arguments{
  \item{d}{A dive (object of class \code{"dive"}).}
  \item{value}{A numeric vector containing elapsed times in minutes.}
}
\details{
  An object of class \code{"dive"} represents a scuba dive. It is
  created by the function \code{\link{dive}}. A dive is defined as a
  series of waypoints occurring at specified depths and times. The
  elapsed time at each waypoint is returned by \code{times.dive}. The
  assignment \code{times.dive(d) <- value} alters these elapsed times,
  provided the new vector \code{value} has the same length as the old one.
}
\value{
  \code{times.dive} returns a numeric vector containing the elapsed
  times at each waypoint, in minutes.
}
\seealso{
  \code{\link{depths.dive}},
  \code{\link{durations.dive}},
  \code{\link{dive}}.
}
\examples{
   d <- dive(c(30,20), c(5,5))
   d
   times.dive(d)
   # stretch time by 10 percent
   times.dive(d) <- 1.1 * times.dive(d)
   d
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
