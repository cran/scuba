\name{Mvalues}
\alias{Mvalues}
\title{
  M-values for Standard Haldane Models
}
\description{
  Contains the surfacing M-values (maximum surfacing nitrogen tensions)
  for several Haldane-type decompression models.
}
\details{
  The surfacing M-value (in the terminology of Workman) is the
  maximum partial pressure of inert gas in a tissue compartment
  which would be tolerated without symptoms of decompression sickness
  at an ambient pressure of 1 atmosphere.

  This dataset gives the surfacing M-values for nitrogen in each tissue,
  in several of the classical decompression models, as
  tabulated by Bookspan (1995, page 23).
}
\format{
  A list with two components called \code{ata} and \code{fsw},
  giving the surfacing M-values in atmospheres absolute
  (\code{ata}) or in equivalent feet of seawater (\code{fsw})
  respectively.
  Each component is a list of three vectors, named
  \tabular{ll}{
    \code{Haldane} \tab J.S. Haldane's original model \cr
    \code{USN} \tab original US Navy model \cr
    \code{DSAT} \tab model for DSAT (PADI) tables
  }
  containing the M-values for the 
  tissues in the relevant model, listed from fastest to slowest.
  The tissues correspond to the tissues in the dataset
  \code{\link{Halftimes}}.
}
\examples{
  data(Mvalues)
  Mvalues$fsw$USN
  Mvalues$ata$DSAT
}
\seealso{
  \code{\link{Halftimes}}
}
\source{
  Bookspan (1995), p 23.
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
}
\keyword{datasets}
