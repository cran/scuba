\name{showstates}
\alias{showstates}
\title{Interactive Display of Diver Saturation}
\description{
  Plots a dive profile and interactively displays the diver's nitrogen
  saturation state at any intermediate stage indicated by the mouse.
}
\usage{
 showstates(d, model="DSAT", relative=TRUE, deco=FALSE, ...)
}
\arguments{
  \item{d}{
    The dive. An object of class \code{"dive"}.
  }
  \item{model}{
    The decompression model. Either an object of class \code{"hm"}
    or a character string identifying a model.
  }
  \item{relative}{
    Logical flag indicating whether to compute relative tissue
    saturations. If \code{FALSE}, tissue saturations are expressed
    as absolute pressures in atmospheres absolute (ata).
    If \code{TRUE} (the default),
    the tissue saturation for each
    compartment is expressed as a fraction of the surfacing M-value
    for the compartment. (Alternatively if \code{deco=TRUE} then tissue
    saturation is expressed as a fraction of the M-value at current depth.)
  }
  \item{deco}{
    Logical flag indicating whether to calculate relative saturations
    for a decompression dive.
    If \code{deco=FALSE}, then relative tissue saturations
    are computed by dividing the absolute tissue saturation by the
    surfacing M-value, as appropriate for a no-decompression dive.
    If \code{deco=TRUE}, then relative tissue saturations are computed
    by dividing the absolute tissue saturation by the M-value at the
    current depth, as appropriate for a decompression dive.
    This argument applies only when \code{relative=TRUE}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.dive}} to control the plot
    of the dive profile.
  }
}
\value{
  The vector of absolute nitrogen tensions (in ata) in the most recently clicked
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

  Note that (by default)
  the bar graph shows the \emph{relative} nitrogen tensions in each
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
    showstates(dive(c(30,20), c(5,3)), "ZH")
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
