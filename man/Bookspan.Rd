\name{Bookspan}
\alias{Bookspan}
\docType{data}
\title{Tissue data from Bookspan's book}
\description{
  Data for several decompression models
  of Haldane type, giving the tissue halftimes for Nitrogen
  and the M-values (maximum tolerated partial pressure of Nitrogen
  for surfacing).
}
\details{
  Bookspan (1995) tabulates the parameters of several different
  decompression models of the classical Haldane type.
  The basic parameters for such a model are
  \describe{
    \item{halftimes:}{the Nitrogen diffusion halftime (in minutes)
      for each tissue in the model}
    \item{M-values:}{the maximum partial pressure of Nitrogen
      in each tissue which would be tolerated without symptoms
      at an ambient pressure of 1 atmosphere (the ``surfacing
      M-value'').
    }
  }
  The data tabulated in Bookspan (1995, page 16) purportedly give the halftimes
  for Haldane's original model, the original US Navy model,
  the DSAT model (basis of PADI repetitive dive planner),
  the OrcaEdge, MicroBrain and Aladin dive computers,
  and the Buehlmann ZH-L12 model.

  The data tabulated in Bookspan (1995, page 23) give the corresponding
  M-values (in fsw) for the US Navy and DSAT models only.

  These parameters are all contained in the dataset \code{Bookspan}.
}
\format{
  The dataset \code{Bookspan} is a list with the following entries:
  \describe{
    \item{Halftimes}{Halftimes data. A list with entries 
      \tabular{ll}{
	\code{Haldane} \tab halftimes for Haldane's original model \cr
	\code{USN} \tab halftimes for US Navy original model \cr
	\code{DSAT} \tab halftimes for DSAT model \cr
	\code{OrcaEdge} \tab halftimes for OrcaEdge computer \cr
	\code{MicroBrain} \tab halftimes for MicroBrain computer \cr
	\code{ZHL12} \tab halftimes for ZH-L12 model
      }
    }
    \item{Mvalues.ata}{surfacing M-values in ata (atmospheres absolute).
      A list with entries
      \tabular{ll}{
	\code{Haldane} \tab M-values for Haldane's original model \cr
	\code{USN} \tab M-values for US Navy original model \cr
	\code{DSAT} \tab M-values for DSAT model
      }
    }
    \item{Mvalues.fsw}{surfacing M-values in fsw (feet of seawater).
      A list with the same structure as \code{Mvalues.fsw}.
    }
  }
}
\source{
   Bookspan (1995), page 16 (halftime data) and page 23 (M-values).
}
\examples{
   data(Bookspan)
   Bookspan$Halftimes$DSAT
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
