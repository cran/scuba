\name{trimix}
\alias{trimix}
\title{Trimix Gas}
\description{
  Create a Trimix gas.
}
\usage{
 trimix(fO2, fHe, fN2)
}
\arguments{
  \item{fO2, fHe, fN2}{
    Fractions (between 0 and 1) of oxygen, helium and nitrogen,
    respectively, in the trimix gas. 
  }
}
\value{
  An object of class \code{"gas"} representing
  Trimix with the specified composition.
}
\details{
  An object of class \code{"gas"} represents a breathing gas.
  Such objects are required for various calculations in the
  \code{scuba} library.
  
  The value returned by \code{trimix} represents Trimix
  (oxygen-nitrogen-helium mixture) with the fractions of oxygen,
  nitrogen and helium specified by the arguments \code{fO2}, \code{fN2},
  \code{fHe} respectively. These fractions should sum to 1.

  If only two of the gas fractions are given, the
  missing fraction will be calculated so that the three values sum to 1.

  There are methods for \code{print} and \code{summary}
  for objects of class \code{"gas"}.
}
\seealso{
  \code{\link{nitrox}},
  \code{\link{air}}
}
\examples{
  # Trimix 18/50
  #            (18\% oxygen, 50\% helium) 
  trimix(0.18, 0.5)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
