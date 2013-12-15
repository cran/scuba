\name{pedro}
\docType{data}
\alias{pedro}
\alias{pedro902}
\alias{pedro903}
\alias{pedro904}
\alias{pedro922}
\alias{pedro943}
\alias{pedro944}
\alias{pedro945}
\alias{pedro946}
\alias{pedro948}
\alias{pedro949}
\alias{pedro950}
\title{
  Real Scuba Dive Profiles
}
\description{
  Dive profiles for 11 different scuba dives.
}
\details{
  The command \code{data(pedro)}
  gives access to 11 datasets named
  \code{pedro902}, 
  \code{pedro903}, 
  \code{pedro904}, 
  \code{pedro922}, 
  \code{pedro943}, 
  \code{pedro944}, 
  \code{pedro945}, 
  \code{pedro946}, 
  \code{pedro948}, 
  \code{pedro949} and
  \code{pedro950}. These are the dive profiles for 11 real scuba dives,
  uploaded from a dive computer.

  Each of the datasets is an object of class \code{"dive"}
  containing the dive profile for one dive. 
  It can be plotted and printed (using \code{\link{plot.dive}}
  and \code{\link{print.dive}}). The nitrogen saturation can be
  computed using \code{\link{haldane}} and the
  cumulative oxygen toxicity using \code{\link{oxtox}}.

  The dives were conducted by Pedro Antonio Neves. Details are
  as follows:
  \tabular{lllll}{
    \emph{Dive} \tab \emph{Max depth} \tab \emph{Gases} \tab \emph{Type}
    \tab \emph{Location}
    \cr
    902 \tab 19 m \tab EAN 32 \tab single \tab Armacao de Pera, Portugal\cr
    903 \tab 27 m \tab EAN 32 \tab single \tab Farilhoes, Portugal\cr
    904 \tab 30 m \tab EAN 32 \tab single \tab Berlenga, Portugal\cr
    922 \tab 15 m \tab EAN 32 \tab multiple \tab Sines, Portugal\cr
    943 \tab 21 m \tab EAN 32 \tab single \tab Armacao de Pera, Portugal\cr
    944 \tab 19 m \tab EAN 32 \tab single \tab Armacao de Pera, Portugal\cr
    945 \tab 22 m \tab EAN 32 \tab single \tab Armacao de Pera, Portugal\cr
    946 \tab 20 m \tab EAN 32 \tab single \tab Armacao de Pera, Portugal\cr
    949 \tab 16 m \tab EAN 32 \tab single \tab Finisterra,Galicia,Spain\cr
    950 \tab 40 m \tab EAN 50, Trimix \tab single \tab Finisterra,Galicia,Spain
  }
  All dives were conducted on Nitrox 32, except \code{pedro950} which was
  conducted using Trimix 18/45 and Nitrox 50.
  All dive profiles consist of a single descent and ascent,
  except for \code{pedro922} which consists of 5 short bounce dives separated
  by short surface intervals. 
}
\format{
  Each of the datasets \code{pedro902}, 
  \code{pedro903}, etc, is an object of class \code{"dive"}.
}
\examples{
   data(pedro)
   summary(pedro950)
   plot(pedro950)
}
\source{
  Pedro Antonio Neves
}
\keyword{datasets}
