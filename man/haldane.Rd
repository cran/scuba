\name{haldane}
\alias{haldane}
\title{Tissue Saturation by Haldane Model}
\description{
  Computes a diver's nitrogen saturation 
  after a dive, as predicted by a Haldane model.
}
\usage{
  haldane(d,
	 halftimes={ data(Halftimes); Halftimes$DSAT}, 
	 prevstate=rep(0.79,length(halftimes)))
}
\arguments{
  \item{d}{
    The dive profile. An object of class \code{\link{dive}}. 
  }
  \item{halftimes}{
    Vector of halftimes for the tissues (compartments) in
    the Haldanean model. Defaults to the values for the
    DSAT (PADI) model.
  }
  \item{prevstate}{
    The diver's nitrogen saturation state at the beginning of the dive,
    for repetitive dives.
    Default is to assume the diver is ``desaturated'' (in equilibrium
    with air at sea level) at start of dive.
  }
}
\value{
  Vector giving the diver's nitrogen saturation state at the end of the
  dive. Each entry is the nitrogen tension (in atmospheres absolute) of
  the compartment.
}
\details{
  Henry's Law is applied 
  to predict the on- and off-gassing of nitrogen
  in each tissue (compartment) of the model.
  The resulting differential equations are solved analytically.
  
  No constraints of any kind are checked. In particular the
  M-values of the model are not used, so it is not guaranteed that the
  model accepts the dive profile as a no-decompression dive.
}
\section{Warnings}{
  Not applicable to altitude dives. Not applicable to
  gas mixtures other than air and nitrox (oxygen-nitrogen mixtures).
}
\seealso{
  \code{\link{dive}},
  \code{\link{oxtox}}
}
\examples{
  # First dive to 25 m for 20 min with safety stop
  d1 <- dive(c(25,20),c(5,5))
  # Evaluate saturation according to DSAT model
  s1 <- haldane(d1)
  # Look at saturation (in ata)
  barplot(s1)
  # Look at saturation (in \% of maximum)
  data(Mvalues)
  barplot(100 * s1/Mvalues$ata$DSAT, ylab="Relative saturation", xlab="Tissues")
  # Evaluate saturation after 2 hour surface interval
  s2 <- haldane(dive(c(0,120)), prevstate=s1)
  # Then after another dive to 18 m for 30 min with safety stop
  s3 <- haldane(dive(c(18, 30),c(5,3)), prevstate=s2)
  # Assess effect of breathing 80\% oxygen at safety stop
  s3o <- haldane(dive(c(18, 30),5, nitrox(0.8), c(5,3)), prevstate=s2)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
