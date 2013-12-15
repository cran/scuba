\name{BuehlmannL16A}
\alias{BuehlmannL16A}
\docType{data}
\title{Decompression model ZH-L16A}
\description{
  The parameters of the decompression model ZH-L16A
  developed by Buehlmann.
}
\details{
  This dataset contains the \emph{tissue parameters} of the
  decompression model ZH-L16A.
  
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

  The dataset \code{BuehlmannL16A} is a table giving the values of
  \eqn{a} and \eqn{b} for each of the 17 tissue compartments.
}
\section{Decompression calculations}{
  The dataset \code{BuehlmannL16A} is just a table of numbers.
  It cannot be used directly for decompression calculations,
  because the \pkg{scuba} package does not recognise it
  as a decompression model (object of class \code{"hm"}).

  To perform decompression calculations with the ZH-L16A model,
  use the functions \code{\link{haldane}} or \code{\link{showstates}}
  with the argument \code{model="ZH-L16A"} (or just \code{model="Z"}).
  
  To convert a table of
  \eqn{a} and \eqn{b} values, such as the table above,
  into an object of class \code{"hm"} representing a decompression model,
  first calculate the M-values \code{M0 = a + 1/b} and \code{dM = 1/b},
  and then pass the values \code{M0}, \code{dM} to the function
  \code{\link{hm}}.
}
\format{
  A data frame giving the following parameters for each compartment:
  \tabular{ll}{
    \code{id} \tab compartment name (1, 1b, 2, \dots, 16) \cr
    \code{halftime} \tab Nitrogen halftime (minutes) \cr
    \code{a} \tab parameter \eqn{a} in \code{ata} (atmospheres absolute) \cr
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
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{datasets}
