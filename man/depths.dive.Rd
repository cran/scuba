\name{depths.dive}
\alias{depths.dive}
\alias{depths.dive<-}
\title{Depths at each waypoint of a dive}
\description{
  Extracts, or alters, the depth at each waypoint during a dive.
}
\usage{
depths.dive(d)
depths.dive(d) <- value
}
\arguments{
  \item{d}{A dive (object of class \code{"dive"}).}
  \item{value}{A numeric vector containing depths in metres.}
}
\details{
  An object of class \code{"dive"} represents a scuba dive profile. It is
  created by the function \code{\link{dive}}. A dive is defined as a
  series of waypoints occurring at specified depths and times. The
  depth at each waypoint is returned by \code{depths.dive}. The
  assignment \code{depths.dive(d) <- value} alters these depths,
  provided the new vector \code{value} has the same length as the old one.
}
\value{
  \code{depths.dive} returns a numeric vector containing the
  depths at each waypoint, in metres.
}
\seealso{
  \code{\link{times.dive}},
  \code{\link{durations.dive}},
  \code{\link{dive}}.
}
\examples{
   d <- dive(c(30,20), c(5,5))
   d
   depths.dive(d)
   # what if we had dived to 35 metres?
   depths.dive(d)[2:3] <- 35
   d
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
