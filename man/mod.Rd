\name{mod}
\alias{mod}
\title{Maximum Operating Depth}
\description{
  Computes the Maximum Operating Depth for a given
  gas mixture.
}
\usage{
 mod(g, ppO2max=1.4)
}
\arguments{
  \item{g}{
    Breathing gas. An object of class \code{"gas"} or a number
    giving the fraction (between 0 and 1) of oxygen in a
    nitrox mixture. 
  }
  \item{ppO2max}{
    maximum permitted partial pressure of oxygen
    in atmospheres absolute
  }
}
\value{
  The maximum operating depth for this gas.
}
\details{
  Computes the maximum depth at which the
  partial pressure of oxygen does not exceed \code{ppO2max}.
}
\section{Warnings}{
  Not applicable to altitude dives. Does not check whether the
  gas would be hypoxic at the surface.
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
  mod(nitrox(0.36))
  # Trimix 15/50
  mod(trimix(0.15, 0.5))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
