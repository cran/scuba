\name{plot.dive}
\alias{plot.dive}
\title{Plot a Dive Profile}
\description{
  Plot a dive profile.
}
\usage{
  \method{plot}{dive}(x, \dots, main=deparse(substitute(x)),
                        key.gases=c("text", "legend", "none"),
                        text.cex=1,
                        text.verticals=TRUE,
                        col.gases=1:length(tanklist(x)),
                        legendpos=c("top","bottom", "left", "right",
                        "topleft", "topright", "bottomleft", "bottomright",
                        "center"))
}
\arguments{
  \item{x}{
    The dive profile. An object of class \code{"dive"}
    created by \code{\link{dive}}. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.default}}.
  }
  \item{main}{
    Overall label for plot. 
  }
  \item{key.gases}{
    String indicating whether the plot should be annotated with labels
    indicating what gas is breathed, and how to do the annotation.
  }
  \item{text.cex}{
    Numeric value controlling the size of text when \code{key.gases="text"}.
  }
  \item{text.verticals}{
    Logical flag indicating whether the plot should
    indicate what gas is breathed during ascent and descent
    between stages, as well as during the main stages of the dive.
  }
  \item{col.gases}{
    Colours to be used to represent each of the breathing gases
    used in the dive. Either a numeric vector of colour codes,
    or a vector of character strings identifying colours.
  }
  \item{legendpos}{
    Character string determining the position of the
    legend, when \code{key.gases="legend"}. Default is \code{"top"}.
  }
}
\value{
  NULL.
}
\details{
  A representation of the dive profile is plotted.
}
\seealso{
  \code{\link{dive}}
}
\examples{
  # Dive to 25 m for 20 min with safety stop
  d1 <- dive(c(25,20),c(5,5))
  # plot it
  plot.dive(d1)
  # Dive to 18 m for 30 min with safety stop, on Nitrox EANx 36
  d3 <- dive(nitrox(0.36), c(18,30), c(5,3))
  plot.dive(d3)
  # Real dive profile
  data(baron)
  d4 <- dive(nitrox(0.30), baron[, 1:2])
  plot(d4, main="Baron Gautsch dive")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
