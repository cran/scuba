\name{hm}
\alias{hm}
\title{Haldane Type Model}
\description{
  Creates a Haldane-type decompression model.
}
\usage{
  hm(HalfT, M0=NULL, dM=NULL, ...,
               N2 = list(HalfT=HalfT, M0=M0, dM=dM),
               He = NULL,
               title="user-defined model",
               cnames=NULL,
               mixrule=NULL) 
}
\arguments{
  \item{HalfT}{
    Vector of nitrogen halftimes (in minutes) for each compartment.
  }
  \item{M0}{
    Optional vector of surfacing M-values (in ata) of nitrogen
    for each compartment.
  }
  \item{dM}{
    Optional vector of gradients for M-values (dimensionless)
    for each nitrogen compartment.
  }
  \item{\dots}{Ignored.}
  \item{N2}{
    An alternative way of specifying all the data for
    Nitrogen. A list with elements labelled
    \code{"HalfT"}, \code{"M0"} and \code{"dM"} giving the
    halftimes, surfacing M-values, and M-value gradients
    for each compartment. 
  }
  \item{He}{
    Data for Helium, if available.
    A list with elements labelled
    \code{"HalfT"}, \code{"M0"} and \code{"dM"} giving the
    Helium halftimes, surfacing M-values, and M-value gradients
    for each compartment.
  }
  \item{title}{
    Optional name of model. A character string.
  }
  \item{cnames}{
    Optional names of compartments. A vector of character strings.
  }
  \item{mixrule}{
    Mixing rule for M-values.
    Either \code{"N2"}, \code{"interpolate"}, or \code{NULL}
    (representing a sensible default). 
  }
}
\value{
  An object of class \code{"hm"}, representing the decompression model.
}
\details{
  This function creates an object of class \code{"hm"},
  which represents a Haldane-type decompression model.
  Objects of this class are needed by the commands \code{\link{ndl}},
  \code{\link{haldane}} and others.

  A Haldane-type decompression model describes the diver's body
  as a set of independent compartments
  connected directly to the breathing gas and governed by classical
  diffusion.

  The argument \code{halfT} specifies the half-times of the 
  nitrogen compartments, in minutes. The length of \code{halfT} implicitly
  determines the number of nitrogen compartments.

  The argument \code{M0}, if present, specifies the surfacing M-values
  for the nitrogen compartments. These are the maximum values of
  nitrogen tissue tension that are tolerated at the surface.
  These values are required in order to plan a no-decompression dive.

  The argument \code{dM}, if present, specifies the rate of increase
  in nitrogen M-values with pressure. The maximum nitrogen tissue
  tension  tolerated at a pressure P atmospheres is \code{M0 + (P-1) * dM}.
  These values are required in order to plan a decompression dive.

  Optionally the model may also allow calculation with Helium diffusion.
  In that case, the argument \code{He} should be a list, with
  components \code{halfT}, \code{M0} and \code{dM}, specifying the
  Helium halftimes, maximum surfacing Helium tensions, and
  Helium gradients, respectively.

  If Helium parameters are included, so that diving with
  \emph{trimix} (Oxygen/Nitrogen/Helium mixture) is permitted,
  then the model must also specify a rule for combining the
  parameters for Helium and Nitrogen to obtain the parameters
  for any trimix gas. This rule is specified by the argument \code{mixrule}.
  Current options are:
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
  The default is \code{mixrule="interpolate"} whenever the Helium
  parameters are specified.
  
  Note that this mixture calculation applies only to the saturation
  parameters \code{M0}, \code{dM} which tell us whether a dive is
  staying within the no-decompression limits. This mixture calculation
  does not affect the Haldane calculations of the gas tensions in the
  diver's body.

  The class \code{"hm"} has methods for \code{print} and
  \code{as.data.frame}.
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
  
  Boycott, A.E. Damant, G.C.C. and Haldane, J.B. (1908)
  The prevention of compressed air illness.
  \emph{Journal of Hygiene} (London) \bold{8}, 342--443.

  Brubakk, A.O. and Neuman, T.S. (eds.) (2003)
  Bennett and Elliott's Physiology and Medicine of Diving.
  5th Edition. Saunders. ISBN 0-7020-2571-2

  Buehlmann, A.A. (1983)
  \emph{Dekompression - Dekompressionskrankheit}.
  Springer-Verlag.
  
  Buehlmann, A.A., Voellm, E.B. and Nussberger, P. (2002)
  \emph{Tauchmedizin}. 5e Auflage. Springer-Verlag.

  Tikvisis, P. and Gerth, W.A. (2003)
  Decompression Theory. 
  In Brubakk and Neuman (2003), Chapter 10.1, pages 419-454.

  Wienke, B.R. (1994)
  \emph{Basic diving physics and applications}.
  Best Publishing Co.
  
  Workman, R.D. (1965) Calculation of decompression schedules
  for nitrogen-oxygen and helium-oxygen dives.
  Research Report 6-65. US Navy Experimental Diving Unit.
  Washington DC. 
}
\seealso{
  \code{\link{pickmodel}} for some standard models of Haldane type;
  \code{\link{ndl}},
  \code{\link{haldane}},
  \code{\link{showstates}}
}
\examples{
  hm(c(10,20,60), rep(2 * 0.79, 3))
  hm(c(10,20,60), rep(2 * 0.79, 3), rep(2 * 0.79, 3))
  hm(c(10,20,60), rep(2 * 0.79, 3), rep(2 * 0.79, 3),
     He=list(HalfT=c(5,10,30), M0=c(1,1,1), dM=c(1,1,1)))
}
\author{\adrian.}
\keyword{utilities}
