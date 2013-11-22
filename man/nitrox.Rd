\name{nitrox}
\alias{nitrox}
\title{Nitrox Mixture}
\description{
  Create a Nitrox gas mixture.
}
\usage{
 nitrox(fO2)
}
\arguments{
  \item{fO2}{
    Fraction (between 0 and 1) of oxygen in the
    nitrox mixture. 
  }
}
\value{
  An object of class \code{"gas"} representing
  Nitrox with the specified fraction of oxygen.
}
\details{
  An object of class \code{"gas"} represents a breathing gas.
  Such objects are required for various calculations in the
  \code{scuba} library.
  
  The value returned by \code{nitrox} represents Nitrox
  (oxygen-nitrogen mixture) with the fraction of oxygen
  specified by the argument \code{fO2}.

  For the standard Nitrox calculations, see \code{\link{ead}}
  and \code{\link{mod}}.
  To determine whether a gas is nitrox, use \code{\link{is.nitrox}}.
  
  There are methods for \code{print} and \code{summary}
  for objects of class \code{"gas"}.
}
\seealso{
  \code{\link{dive}},
  \code{\link{ead}},
  \code{\link{mod}},
  \code{\link{air}},
  \code{\link{is.nitrox}}
}
\examples{
  # Nitrox I (32\% oxygen) 
  nitrox(0.32)
  # Nitrox II (36\% oxygen)
  nitrox(0.36)
  #
  g <- nitrox(0.50)
  is.nitrox(g)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
