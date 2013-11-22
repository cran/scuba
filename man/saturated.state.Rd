\name{saturated.state}
\alias{saturated.state}
\title{Saturated Tissue State}
\description{
  Computes the amount of inert gas in each compartment of a diver
  who has reached equilibrium with the breathing gas
  at the specified depth.
  The default is a fresh diver with no diving history.
}
\usage{
saturated.state(model="D", depth = 0, g = air)
}

\arguments{
  \item{model}{
    The decompression model. An object of class \code{"hm"}.
    Defaults to the DSAT (PADI) model.
  }
  \item{depth}{
    The diving depth at which saturation is achieved.
    A single number.
  }
  \item{g}{
    The breathing gas. An object of class \code{"gas"}.
  }
}
\details{
  A Haldane-type decompression model describes the diver's body
  as a set of independent compartments
  connected directly to the breathing gas and governed by classical
  diffusion.

  This command computes the amount of inert gas in each compartment of a diver
  who has reached equilibrium with the breathing gas
  at the specified depth. 

  The default is a fresh diver (in equilibrium with air at the surface
  at sea level).
}
\value{
  A data frame, with one row for each compartment of the model,
  and one column for each inert gas in the model (\code{N2} and/or
  \code{He}). Entries in the data frame are tissue tensions in ata.
}
\seealso{
  \code{\link{hm}},
  \code{\link{haldane}},
  \code{\link{nitrox}},
  \code{\link{trimix}}
}

\examples{
  fresh <- saturated.state()
  deep  <- saturated.state("ZH", 80, trimix(0.15, 0.5))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
