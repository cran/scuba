\name{Workman65}
\alias{Workman65}
\title{Decompression model of Workman 1965}
\description{
  The parameters of the decompression model
  due to Workman, 1965.
}
\details{
  The decompression model developed for the US Navy by Workman (1965)
  is a pure diffusion (Haldane-type) model consisting of 9 tissue compartments.
  Nitrogen in the breathing gas diffuses in and out of each compartment
  at a rate governed by the halftime for that compartment.
  The maximum tolerable nitrogen tension in a compartment,
  when the diver is at depth \eqn{D} msw (metres of sea water), is
  \eqn{M = M_0 + A \times D}{M = M0 + A * D}
  msw, where \eqn{M_0}{M0} and \eqn{A} are values intrinsic to the tissue.
  The value \eqn{M_0}{M0} is called the `surfacing M-value'
  since it is the maximum nitrogen tension that is allowed
  in each compartment during a no-decompression dive at sea level.
}
\format{
  A data frame giving the following parameters for each compartment:
  \tabular{ll}{
    \code{halftime} \tab Nitrogen halftime (minutes) \cr
    \code{M0} \tab surfacing M-value parameter \eqn{M_0}{M0} in msw (metres of seawater) \cr
    \code{A} \tab slope parameter \eqn{A} (dimensionless)
  }
}
\source{
  Tikvisis and Gerth (2003), Table 10.1.1, page 440.
}
\references{
  Brubakk, A.O. and Neuman, T.S. (eds.) (2003)
  Bennett and Elliott's Physiology and Medicine of Diving.
  5th Edition. Saunders. ISBN 0-7020-2571-2

  Tikvisis, P. and Gerth, W.A. (2003)
  Decompression Theory. 
  In Brubakk and Neuman (2003), Chapter 10.1, pages 419-454.

  Workman, R.D. (1965) Calculation of decompression schedules
  for nitrogen-oxygen and helium-oxygen dives.
  Research Report 6-65. US Navy Experimental Diving Unit.
  Washington DC. 
}
\seealso{
  \code{\link{BuehlmannL16A}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{datasets}
