\name{ascent}
\alias{ascent}
\title{Ascent Rate or Time}
\description{
  Specify an Ascent Rate or Ascent Time.
}
\usage{
 ascent(speed=NULL, time=NULL)
}
\arguments{
  \item{speed}{
    Ascent rate in metres per minute.
    Incompatible with \code{time}.
  }
  \item{time}{
    Total ascent time in minutes.
    Incompatible with \code{speed}.
  }
}
\value{
  An object of class \code{"rate"} representing
  the ascent rate or ascent time.
}
\details{
  An object of class \code{"rate"} represents a diver's
  rate of ascent or descent, or the total time taken to
  ascend or descend. Such objects are useful in describing
  dives, especially in the function \code{\link{dive}}.

  The value returned by \code{ascent} represents an ascent rate
  or ascent time. Exactly one of the arguments \code{speed} and
  \code{time} should be given. If \code{speed} is given, this specifies
  a fixed rate of ascent, in metres per minute. If \code{time} is given,
  this specifies a fixed total time of ascent, in minutes.
}
\seealso{
  \code{\link{dive}},
  \code{\link{descent}}
}
\examples{
  # Ascend at 18 metres/minute
  ascent(18)
  # Ascend in exactly 3 minutes
  ascent(time=3)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
