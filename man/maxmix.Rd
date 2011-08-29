\name{maxmix}
\alias{maxmix}
\title{Optimal Nitrox Mixture For Given Depth}
\description{
  Computes the optimal nitrox mixture for a
  given maximum depth. 
}
\usage{
 maxmix(depth, ppO2max=1.4)
}
\arguments{
  \item{depth}{
    the maximum depth, in metres
  }
  \item{ppO2max}{
    maximum permitted partial pressure of oxygen
    in atmospheres absolute
  }
}
\value{
  The optimal nitrox mixture. An object of class \code{"gas"}.
}
\details{
  Computes the maximum fraction of oxygen in a nitrox mixture
  subject to the constraint that the partial pressure of oxygen
  does not exceed \code{ppO2max} atmospheres. 
}
\section{Warnings}{
  Not applicable to altitude dives. Not applicable to
  gas mixtures other than nitrox (oxygen-nitrogen mixtures).
}
\seealso{
  \code{\link{mod}},
  \code{\link{ead}},
  \code{\link{eadtable}}
}
\examples{
  # 30 metres
  maxmix(30)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
