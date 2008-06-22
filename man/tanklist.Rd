\name{tanklist}
\alias{tanklist}
\alias{tanklist<-}
\title{Extract or Change the Breathing Gas Tanks in a Dive}
\description{
  Extracts or modifies the list of breathing gases
  in a dive object.
}
\usage{
tanklist(d)
tanklist(d) <- value
}
\arguments{
  \item{d}{The dive (an object of class \code{"dive"}).}
  \item{value}{The new tank list. A list, whose entries are gases
    (objects of class \code{"gas"})
  }
}
\details{
  An object of class \code{"dive"} represents a scuba dive,
  including information about depth, time and gas breathed at each
  stage of the dive. These objects are created by the function
  \code{\link{dive}}.

  The tank list of a dive object is a list of the breathing gases
  that were used during the dive. The function \code{tanklist} returns
  this list.

  If a new value is assigned to the tank list of a dive \code{d},
  then \code{d} is changed. The new dive \code{d} is conducted
  to the same depths and times as the old \code{d}, but with different
  gases. 
}
\value{
  The value of \code{tanklist(d)} is a list whose elements are
  gases (objects of class \code{"gas"}). 
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{dive}},
  \code{\link{nitrox}}
}
\examples{
   d <- dive(air, c(30,4), 5, nitrox(0.5), c(5,10))
   d
   tanklist(d)
   tanklist(d) <- list(air, nitrox(0.8))
   d
   tanklist(d) <- list(travel=air, deco=nitrox(0.8))
   d
}
\keyword{utilities}
