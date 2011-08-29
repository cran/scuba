\name{showstates}
\alias{showstates}
\title{Interactive Display of Diver Saturation}
\description{
  Plots a dive profile and interactively displays the diver's nitrogen
  saturation state at any intermediate stage indicated by the mouse.
}
\usage{
 showstates(d, model="DSAT")
}
\arguments{
  \item{d}{
    The dive. An object of class \code{"dive"}.
  }
  \item{model}{
    The decompression model. Either an object of class \code{"hm"}
    or a character string identifying a model.
  }
}
\value{
  The vector of nitrogen tensions (in ata) in the most recently clicked
  stage.
}
\details{
  An object of class \code{"dive"} represents a scuba dive.
  It is created by the command \code{\link{dive}}.

  This function plots the dive using \code{\link{plot.dive}}
  then waits for the user to click on the plot. This click selects
  a time during the dive. The algorithm computes the nitrogen tensions
  in the tissues at the indicated time, using \code{\link{haldane}},
  and plots them as a bar graph. The cumulative oxygen toxicity is also
  computed using \code{\link{oxtox}}.

  Note that the bar graph shows the \emph{relative} nitrogen tensions in each
  compartment, that is, the tissue nitrogen tension divided by the
  "surfacing M-value" (tissue maximum nitrogen tension for a
  no-decompression dive at sea level). 

  The argument \code{model} determines the decompression model.
  It should be either an object of class \code{"hm"} (created by
  \code{\link{hm}}) or a character string matching one of the
  options in \code{\link{pickmodel}}.
}
\seealso{
  \code{\link{dive}},
  \code{\link{haldane}},
  \code{\link{oxtox}}
}
\examples{
  \dontrun{
    showstates(dive(c(30,20)), "ZH")
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
