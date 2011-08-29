\name{is.nitrox}
\alias{is.nitrox}
\title{Recognise a Nitrox Gas}
\description{
  Determines whether a given gas is nitrox.
}
\usage{
  is.nitrox(g)
}
\arguments{
  \item{g}{
    An object of class \code{"gas"}.
  }
}
\details{
  An object of class \code{"gas"} represents a breathing gas.
  Such objects are required for various calculations in the
  \code{scuba} library.

  The function \code{is.nitrox} expects its argument \code{g}
  to be a gas object. It returns \code{TRUE} if 
  \code{g} is a nitrox gas (a mixture of oxygen and nitrogen),
  and \code{FALSE} otherwise.
}
\seealso{
  \code{\link{nitrox}}
}
\examples{
  v <- nitrox(0.50)
  is.nitrox(v)
  is.nitrox(air)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
