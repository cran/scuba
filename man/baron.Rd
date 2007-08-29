\name{baron}
\alias{baron}
\title{
  Real Scuba Dive Profile 
}
\description{
  Contains the dive profile data for a scuba dive
  on the \emph{Baron Gautsch} wreck, uploaded from
  a dive computer. 
}
\details{
  This dataset contains the dive profile data for a real scuba dive
  uploaded from a dive computer. It gives the
  depth, water temperature, breathing gas cylinder pressure,
  calculated residual bottom time, and calculated workload,
  recorded every 4 seconds throughout the dive.

  The dive site was the wreck of the vessel \emph{Baron Gautsch} in Croatia.
  The dive, to a maximum depth of 40 metres,
  with total dive time of 40 minutes,
  was conducted on nitrox (EAN 30) by Vittorio Broglio.

  The examples show how to convert the depth-time data to a \code{\link{dive}}
  object.
}
\format{
  A data frame with the following columns:
  \tabular{lll}{
    \code{time} \tab character \tab elapsed time (minutes:seconds) \cr
    \code{depth} \tab numeric \tab depth (metres) \cr
    \code{temp} \tab numeric \tab temperature (C) \cr
    \code{bar} \tab numeric \tab breathing gas pressure (bar) \cr
    \code{RBT} \tab integer \tab residual bottom time (min) \cr
    \code{WL} \tab integer \tab workload 
  }
}
\examples{
  data(baron)
  b <- dive(nitrox(0.30), baron[,1:2])
}
\source{
  Vittorio Broglio
}
\keyword{datasets}
