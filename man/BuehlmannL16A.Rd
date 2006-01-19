\name{BuehlmannL16A}
\alias{BuehlmannL16A}
\title{Decompression model ZH-L16A}
\description{
  The parameters of the decompression model ZH-L16A
  developed by Buehlmann.
}
\details{
  Buehlmann's decompression model ZH-L16A is a pure
  diffusion (Haldane-type) model consisting of 17 tissue compartments.
  Nitrogen in the breathing gas diffuses in and out of each compartment
  at a rate governed by the halftime for that compartment.
  The maximum tolerable nitrogen tension in a compartment,
  when the ambient pressure is \eqn{P} atmospheres absolute,
  is
  \deqn{P_{\mbox{tol}} = a + \frac{P}{b}}{Ptol = a + P/b}
  where \eqn{a} and \eqn{b} are values intrinsic to the tissue.
  Thus, for a no-decompression dive at sea level, the nitrogen tension
  in each compartment must not exceed the value \eqn{a + 1/b} for
  that compartment.
}
\format{
  A data frame giving the following parameters for each compartment:
  \tabular{ll}{
    \code{id} \tab compartment name (1, 1b, 2, \dots, 16) \cr
    \code{halftime} \tab Nitrogen halftime (minutes) \cr
    \code{a} \tab parameter \eqn{a} in ata (atmospheres absolute) \cr
    \code{b} \tab parameter \eqn{b} (dimensionless)
  }
}
\source{
  Buehlmann et al (2003), Table 25, page 158.
}
\references{
  Buehlmann, A.B., Voellm, E.B. and Nussberger, P. (2002)
  Tauchmedizin. 5th edition. Springer-Verlag, Berlin.
  ISBN 3-540-42979-4.
}
\seealso{
  \code{\link{Workman65}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{datasets}
