\name{ppO2}
\alias{ppO2}
\title{Oxygen Partial Pressures}
\description{
  Computes the partial pressure of oxygen at each stage during a dive.
}
\usage{
  ppO2(d)
}
\arguments{
  \item{d}{
    The dive profile. An object of class \code{"dive"}.
  }
}
\value{
  A data frame with columns \code{time} and \code{ppO2}.
}
\details{
  This function computes the partial pressure of oxygen at each stage of the
  dive profile. The result is a data frame with two columns: \code{time}
  giving the elapsed time, and \code{ppO2} giving the partial pressure
  of oxygen (in atmospheres) in the current breathing gas at the current depth.
}
\seealso{
  \code{\link{oxtox}},  \code{\link{ead}},
  \code{\link{mod}},
  \code{\link{maxmix}}
}
\examples{
  data(pedro)
  plot(ppO2(pedro902), type="l")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
