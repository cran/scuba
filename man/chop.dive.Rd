\name{chop.dive}
\alias{chop.dive}
\title{Extract Part of a Dive Profile}
\description{
  Extracts only part of a dive profile.
}
\usage{
chop.dive(d, t0=0, t1=max(times.dive(d)))
}
\arguments{
  \item{d}{The dive. An object of class \code{"dive"}.}
  \item{t0, t1}{Elapsed times (in minutes) of the
    start and end of the period which should be extracted.}
}
\details{
  The argument \code{d} should be an object of class \code{"dive"}
  representing a dive profile.

  This command extracts the part of the dive
  profile that starts at time \code{t0} minutes and ends at
  time \code{t1} minutes, and returns it as
  an object of class \code{"dive"}. The clock is adjusted so that the
  new dive profile starts at time \code{0} and ends at time \code{t1-t0}.

  Note that the resulting dive profile does not start and end at the surface:
  the dive profile is
  simply chopped off at elapsed time \code{t1}. If you want the dive
  to start and end at the surface, execute \code{dive(chop.dive(d, t0, t1))}.
}
\value{
  An object of class \code{"dive"} starting 
  at time \code{0} and ending at time \code{t1-t0}. 
}
\seealso{
  \code{\link{dive}}
}
\examples{
   d <- dive(c(30,20), c(5,5))
   # fragment of dive up to 10 minutes
   chop.dive(d, 0, 10)
   # dive aborted at 10 minutes
   dive(chop.dive(d, 0, 10))
}
\author{\adrian.}
\keyword{utilities}
