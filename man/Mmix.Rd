\name{Mmix}
\alias{Mmix}
\alias{M0mix}
\title{
  Compute M-values for a Mixture of Inert Gases
}
\description{
  Computes the maximum tolerated inert gas tension at depth
  (M-value) or at the surface (surfacing M-value, \code{M0})
  for a given gas.
}
\usage{
  M0mix(model, fN2, fHe, mixrule=NULL)

  Mmix(model, depth, fN2, fHe, mixrule=NULL)
}
\arguments{
  \item{model}{
    The decompression model.
    Either an object of class \code{"hm"} (see \code{\link{hm}})
    or a character string giving the name of one of the
    existing models (see \code{\link{pickmodel}}).
  }
  \item{depth}{
    Numeric value or numeric vector giving the depth or depths
    in metres of seawater.
  }
  \item{fN2}{
    Fraction of nitrogen in the breathing gas.
    Number between 0 and 1.
    Alternatively this can be a \code{gas} object,
    and the values of \code{fN2,fHe} will be extracted from it.
  }
  \item{fHe}{
    Fraction of helium in the breathing gas.
    Number between 0 and 1.
  }
  \item{mixrule}{
    Optional string specifying the mixture rule:
    either \code{"N2"} or \code{"interpolate"}.
    See Details.
    If this argument is absent, the mixture rule
    is taken from the \code{model}. If this argument is present,
    it overrides the mixture rule in \code{model}.
  }
}
\details{
  \code{M0mix} computes the maximum tolerated inert gas tension,
  in each tissue compartment, for a diver returning to the surface
  (known as the surfacing M-value \code{M0}). This is used for
  planning no-decompression dives.

  \code{Mmix} computes the maximum tolerated inert gas tension,
  in each tissue compartment, for a diver at a specified depth,
  called the M-value. This is used for planning decompression dives.

  A Haldane-type decompression model only specifies
  these M-values for a single inert gas (i.e. only nitrogen
  or only helium). When the breathing gas is trimix,
  the inert gas is a mixture of nitrogen and helium, and 
  the M-values for nitrogen and helium must somehow be combined
  to obtain M-values relevant to the mixed gas. 

  The rule for combining the M-values is specified
  by the argument \code{mixrule}. Current options are:
  \describe{
    \item{\code{"N2"}}{
      Ignore the Helium parameters; pretend that Helium is Nitrogen.
      Combine Nitrogen and Helium into a single inert gas,
      and take the parameters \code{M0}, \code{dM} for this gas
      to be the parameters \code{M0}, \code{dM} for Nitrogen.
    }
    \item{\code{"interpolate"}}{
      Apply Buehlmann's (1983, 2002) interpolation rule. Convert the
      parameters \code{M0}, \code{dM} to the Buehlmann parameters
      \code{a = M0 - dM} and \code{b = 1/dM}. For a mixture of
      Nitrogen and Helium, calculate the \code{a,b} values by
      linear interpolation between the values for Nitrogen and Helium
      according to the gas fractions. Then convert from \code{a,b}
      back to \code{M0, dM}.
    }
  }
  If the argument \code{mixrule} is missing or \code{NULL}, then the
  default mixture rule is taken from the \code{model}
  (since every Haldane model has a mixture rule,
  as explained in \code{\link{hm}}).
}
\value{
  A matrix, with one column for each compartment in the model,
  and one row for each entry of \code{fN2} (and \code{fHe} and
  \code{depth}).
  Entries in the matrix are M-values in atmospheres absolute (ata).
}
\author{\adrian.}
\seealso{
  \code{\link{hm}}, 
  \code{\link{pickmodel}} for specifying models.

  \code{\link{haldane}} for computing tissue saturations.
}
\examples{
  # Trimix 18/50, surfacing M-values
  M0mix("Z", trimix(0.18, 0.5))

  # Trimix 18/50, M-values at 0, 10, 20 metres
  Mmix("Z", c(0,10,20), trimix(0.18, 0.5))
}
\keyword{utilities}
