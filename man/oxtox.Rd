\name{oxtox}
\alias{oxtox}
\title{Pulmonary Oxygen Toxicity}
\description{
  Computes pulmonary oxygen toxicity dose for a given
  dive profile and breathing gas (air or nitrox).
}
\usage{
  oxtox(d, progressive=FALSE)
}
\arguments{
  \item{d}{
    The dive profile. An object of class \code{"dive"}.
  }
  \item{progressive}{
    Logical flag. If \code{FALSE}, the total oxygen toxicity from the
    dive is calculated. If \code{TRUE}, the cumulative oxygen toxicity
    at each time point during the dive is calculated.
  }
}
\value{
  If \code{progressive=FALSE}, a single numeric value,
  the total pulmonary oxygen toxicity dose in OTU.

  If \code{progressive=TRUE}, a numeric vector containing the 
  cumulative pulmonary oxygen toxicity dose in OTU at each
  time point during the dive. The corresponding time points can be
  extracted by \code{\link{times.dive}}.
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
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.

  Brubakk, A.O. and Neuman, T.S. (eds.) (2003)
  Bennett and Elliott's Physiology and Medicine of Diving.
  5th Edition. Saunders. ISBN 0-7020-2571-2

  Buehlmann, A.A., Voellm, E.B. and Nussberger, P. (2002)
  \emph{Tauchmedizin}. 5e Auflage. Springer-Verlag.

  Wienke, B.R. (1994)
  \emph{Basic diving physics and applications}.
  Best Publishing Co.
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
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
