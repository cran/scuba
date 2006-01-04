\name{showstates}
\alias{showstates}
\title{Interactive Display of Diver Saturation}
\description{
  Plots a dive profile and interactively displays the diver's nitrogen
  saturation state at any intermediate stage indicated by the mouse.
}
\usage{
 showstates(d, tissues="DSAT")
}
\arguments{
  \item{d}{
    The dive. An object of class \code{"dive"}.
  }
  \item{tissues}{
    Either \code{"DSAT"} or \code{"USN"} indicating the model.
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
}
\seealso{
  \code{\link{dive}},
  \code{\link{haldane}},
  \code{\link{oxtox}}
}
\examples{
  \dontrun{
    showstates(dive(c(30,20)))
  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
