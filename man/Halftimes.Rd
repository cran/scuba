\name{Halftimes}
\alias{Halftimes}
\title{
  Half-times for Standard Haldane Models
}
\description{
  Contains the `tissue' (compartment) halftimes for several Haldane-type
  decompression models.
} 
\format{
  A list of vectors. The components are named
  \tabular{ll}{
    \code{Haldane} \tab J.S. Haldane's original model halftimes\cr
    \code{USN} \tab US Navy \cr
    \code{DSAT} \tab model for DSAT (PADI) tables \cr
    \code{OrcaEdge} \tab model for Dacor Orca Edge dive computer \cr
    \code{MicroBrain} \tab model for Beuchat MicroBrain dive computer \cr
    \code{Aladin} \tab model for Uwatec Aladin dive computer \cr
    \code{ZHL12} \tab Buehlmann's ZH-L12 model 
  }
  Each vector gives the nitrogen halftimes in minutes
  for the tissues in the relevant model, listed from fastest to slowest.
}
\examples{
  data(Halftimes)
  Halftimes$Haldane
  Halftimes[["USN"]]
}
\source{
  Bookspan (1995) p 16.
}
\seealso{
  \code{\link{Mvalues}}
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
}
\keyword{datasets}

