\name{descent}
\alias{descent}
\title{Descent Rate or Time}
\description{
  Specify an Descent Rate or Descent Time.
}
\usage{
 descent(speed=NULL, time=NULL)
}
\arguments{
  \item{speed}{
    Descent rate in metres per minute.
    Incompatible with \code{time}.
  }
  \item{time}{
    Total descent time in minutes.
    Incompatible with \code{speed}.
  }
}
\value{
  An object of class \code{"rate"} representing
  the descent rate or descent time.
}
\details{
  An object of class \code{"rate"} represents a diver's
  rate of ascent or descent, or the total time taken to
  ascend or descend. Such objects are useful in describing
  dives, especially in the function \code{\link{dive}}.

  The value returned by \code{descent} represents a descent rate
  or descent time. Exactly one of the arguments \code{speed} and
  \code{time} should be given. If \code{speed} is given, this specifies
  a fixed rate of descent, in metres per minute. If \code{time} is given,
  this specifies a fixed total time of descent, in minutes.
}
\seealso{
  \code{\link{dive}},
  \code{\link{descent}}
}
\examples{
  # Descend at 30 metres/minute
  descent(30)
  # Descend in exactly 3 minutes
  descent(time=3)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
