\name{mod}
\alias{mod}
\title{Maximum Operating Depth}
\description{
  Computes the Maximum Operating Depth for a given
  nitrox mixture.
}
\usage{
 mod(g, ppO2max=1.4)
}
\arguments{
  \item{g}{
    Breathing gas. An object of class \code{"gas"} or a number
    giving the fraction (between 0 and 1) of oxygen in the
    nitrox mixture. 
  }
  \item{ppO2max}{
    maximum permitted partial pressure of oxygen
    in atmospheres absolute
  }
}
\value{
  The maximum operating depth for this nitrox mixture.
}
\details{
  Computes the maximum depth at which the
  partial pressure of oxygen does not exceed \code{ppO2max}.
}
\section{Warnings}{
  Not applicable to altitude dives. Not applicable to
  gas mixtures other than nitrox (oxygen-nitrogen mixtures).
}
\seealso{
  \code{\link{ead}},
  \code{\link{eadtable}},
  \code{\link{maxmix}}
}
\examples{
  # Nitrox I (32\% oxygen) 
  mod(0.32)
  # Nitrox II (36\% oxygen)
  mod(0.36)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
