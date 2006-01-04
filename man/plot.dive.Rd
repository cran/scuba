\name{plot.dive}
\alias{plot.dive}
\title{Plot a Dive Profile}
\description{
  Plot a dive profile.
}
\usage{
  plot.dive(x, \dots, main, verticals=TRUE)
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
  \item{verticals}{
    Flag indicating whether the plot should
    indicate what gas is breathed during ascent and descent
    between stages, as well as during the main stages of the dive.
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
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
