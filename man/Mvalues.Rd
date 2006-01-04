\name{Mvalues}
\alias{Mvalues}
\alias{Mvalues.ata}
\alias{Mvalues.fsw}
\title{
  M-values for Standard Haldane Models
}
\description{
  Contains the surface M-values (maximum surfacing nitrogen tensions)
  for several Haldane-type decompression models.
} 
\format{
  A list of vectors. The components are named
  \describe{
    \item{USN}{US Navy}
    \item{DSAT}{model for DSAT (PADI) tables}
  }
  Each vector gives the M-values in atmospheres absolute
  (\code{Mvalues.ata}) or in equivalent feet of seawater (\code{Mvalues.fsw})
  for the tissues in the relevant model, listed from fastest to slowest.
}
\examples{
  Mvalues.ata[["USN"]]
  Mvalues.ata[["DSAT"]]
  Mvalues.fsw[["USN"]]
  Mvalues.fsw[["DSAT"]]
}
\source{Jolie Bookspan's book}
\section{Notes}{
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
}
\keyword{datasets}
