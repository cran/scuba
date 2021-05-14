\name{END}
\alias{END}
\title{Equivalent Air Depth}
\description{
  Computes the Equivalent Narcotic Depth for a Trimix mixture at
  given depth.
}
\usage{
 END(depth, g)
}
\arguments{
  \item{depth}{
    depth of dive, in metres.
    A numeric value or a vector of values.
  }
  \item{g}{
    The breathing gas. An object of class \code{"gas"}.
  }
}
\value{
  Numeric value or vector giving the Equivalent Narcotic Depth
  in metres for each value of \code{depth}.
}
\details{
  Applies the standard formula for equivalent narcotic depth:
  the depth at which compressed air would have the same
  partial pressure of narcotic gases as the specified gas \code{g}
  at the specified \code{depth}
  (assuming that both Nitrogen and Oxygen are narcotic). 
}
\section{Warnings}{
  This function does not check whether the breathing gas would be safe
  (it could be hypoxic or toxic at the depth in question).
  It is not applicable to altitude dives. 
}
\seealso{
  \code{\link{ead}}
}
\examples{
  # Trimix 18/50
  #            (18\% oxygen, 50\% helium) 
  END(30, trimix(0.18, 0.5))
}
\author{\adrian.}
\keyword{utilities}
