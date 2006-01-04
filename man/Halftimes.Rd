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
  \describe{
    \item{Haldane}{J.S. Haldane's original model halftimes}
    \item{USN}{US Navy}
    \item{DSAT}{model for DSAT (PADI) tables}
    \item{OrcaEdge}{model for Dacor Orca Edge dive computer}
    \item{MicroBrain}{model for Beuchat MicroBrain dive computer}
    \item{Aladin}{model for Uwatec Aladin dive computer}
    \item{Buehlmann.ZHL.12}{Buehlmann's ZHL.12 model}
  }
  Each vector gives the halftimes in minutes
  for the tissues in the relevant model, listed from fastest to slowest.
}
\examples{
  Halftimes[["Haldane"]]
  Halftimes[["USN"]]
  Halftimes[["DSAT"]]
  Halftimes[["OrcaEdge"]]
  Halftimes[["MicroBrain"]]
  Halftimes[["Aladin"]]
  Halftimes[["Buehlmann.ZHL.12"]]
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

