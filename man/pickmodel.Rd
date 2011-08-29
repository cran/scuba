\name{pickmodel}
\alias{pickmodel}
\title{Standard Decompression Models} 
\description{
  Selects one of a list of standard Haldane-type
  decompression models, identified by the name of the model.
}
\usage{
pickmodel(model)
}
\arguments{
  \item{model}{The name of the model. A character string.}
}
\details{
  A Haldane-type decompression model describes the diver's body
  as a set of independent compartments connected directly to the
  breathing gas and governed by classical diffusion. 

  This function is a convenient way to select one of the
  commonly used decompression models of Haldane type.

  The argument \code{model} is partially matched to
  one of the following:
  \describe{
    \item{"Haldane"}{
      The original Haldane (1908) model
      with 5 compartments.
    }
    \item{"USN"}{
      The first model used by the US Navy,
      with 6 compartments.
    }
    \item{"DSAT"}{
      The DSAT model, with 8 compartments,
      underlying the PADI recreational
      dive planner tables.
    }
    \item{"Workman65"}{
      The model of Workman (1965) with 9 compartments.
    }
    \item{"ZH-L16A"}{
      The `theoretical' model ZH-L16A of Buehlmann
      with 17 compartments.
    }
  }
  The result is an object (of class \code{"hm"})
  in the format required for use by commands such as
  \code{\link{haldane}}, \code{\link{ndl}} and \code{\link{showstates}}.

  The data for the USN and Haldane models are taken from
  Wienke (1994), page 127 ff. The DSAT model is
  taken from   Bookspan (1995), pp 16 and 23.
  The Workman 65 model is taken from
  Tikvisis and Gerth (2003), Table 10.1.1, page 440.
  The ZH-L16A parameters
  were taken from Buehlmann et al (2002) Table 25, page 158.
}
\value{
  An object of class \code{"hm"} representing a Haldane-type
  decompression model.
}
\section{Note}{
  If the argument \code{model} is missing, then the list of available
  options for \code{model} is printed.
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
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
\seealso{
  \code{\link{ndl}},
  \code{\link{haldane}},
  \code{\link{showstates}}
}
\examples{
   pickmodel("Z")
   ndl(30, model=pickmodel("Z"))
}
\keyword{utilities}
