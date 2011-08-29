\name{print.dive}
\alias{print.dive}
\title{Print a Dive Profile}
\description{
  Print a dive profile.
}
\usage{
  \method{print}{dive}(x, \dots, seconds=TRUE)
}
\arguments{
  \item{x}{
    The dive profile. An object of class \code{"dive"}
    created by \code{\link{dive}}. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{print.default}}.
  }
  \item{seconds}{
    Flag indicating whether to print elapsed time to the
    nearest second (if \code{TRUE}) or to the nearest minute (if \code{FALSE}).
  }
}
\value{
  NULL.
}
\details{
  A representation of the dive profile is printed.
}
\seealso{
  \code{\link{dive}},
  \code{\link{plot.dive}}
}
\examples{
  # Dive to 25 m for 20 min with safety stop
  dive(c(25,20),c(5,5))
  # Dive to 18 m for 30 min with safety stop, on Nitrox EANx 36
  dive(nitrox(0.36), c(18,30), c(5,3))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
