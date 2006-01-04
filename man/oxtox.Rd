\name{oxtox}
\alias{oxtox}
\title{Pulmonary Oxygen Toxicity}
\description{
  Computes pulmonary oxygen toxicity dose for a given
  dive profile and breathing gas (air or nitrox).
}
\usage{
  oxtox(d)
}
\arguments{
  \item{d}{
    The dive profile. An object of class \code{"dive"}.
  }
}
\value{
  The pulmonary oxygen toxicity dose in OTU.
}
\details{
  Computes the total dose of pulmonary oxygen toxicity
  from the given dive profile, by
  \deqn{
    \int_0^T \left( \frac{\mbox{ppO}_2 - 0.5}{0.5} \right)^{0.83} \; dt
  }{
    integral of ((ppO2 - 0.5)/0.5)^0.83  dt
  }
  The maximum tolerable dose per day is usually reckoned as 1500
  OTU. Allowing 650 OTU for recompression therapy implies a working
  maximum of 850 OTU per day. 
}
\section{Warnings}{
  Not applicable to altitude dives. Not applicable to
  gas mixtures other than nitrox (oxygen-nitrogen mixtures).
}
\seealso{
  \code{\link{ead}},  \code{\link{eadtable}},
  \code{\link{mod}},
  \code{\link{maxmix}}
}
\examples{
  # Nitrox II (36\% oxygen) at 30 metres for 27 minutes
  d <- dive(nitrox(0.36), c(30,27))
  oxtox(d)

  # Same as above, followed by safety stop on 100\% oxygen 
  d <- dive(nitrox(0.36), c(30,27),5, nitrox(1), c(5,5))
  oxtox(d)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
