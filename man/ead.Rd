\name{ead}
\alias{ead}
\alias{eadtable}
\title{Equivalent Air Depth}
\description{
  Computes the Equivalent Air Depth for a Nitrox mixture at
  given depths.
}
\usage{
 ead(depth, g)
 eadtable(g, ppO2max=1.4)
}
\arguments{
  \item{depth}{
    depth of dive, in metres.
  }
  \item{g}{
    The breathing gas. An object of class \code{"gas"}
    or a number specifying the fraction (between 0 and 1) of oxygen in the
    nitrox mixture. 
  }
  \item{ppO2max}{
    maximum permitted partial pressure of oxygen
    in ata.
  }
}
\value{
  For \code{ead}, the Equivalent Air Depth in metres for each
  value of \code{depth};
  for \code{eadtable}, a data frame containing depths and Equivalent
  Air Depth values, for a range of depths up to the maximum permitted
  by \code{ppO2max}.
}
\details{
  Applies the standard formula for equivalent air depth.
}
\section{Warnings}{
  This function does not check whether the breathing gas would be safe
  (it could be hypoxic or toxic at the depth in question).

  \code{ead} returns negative values for sufficiently shallow depths.
  In \code{eadtable} these values are blank.
  
  Not applicable to altitude dives. Not applicable to
  gas mixtures other than nitrox (oxygen-nitrogen mixtures).
}
\seealso{
  \code{\link{mod}},
  \code{\link{maxmix}}
}
\examples{
  # Nitrox I (32\% oxygen) at 20 metres
  ead(20, 0.32)
  # Nitrox I table of EAD's
  eadtable(0.32)
  # Nitrox II (36\% oxygen) at a range of depths
  ead(10:25, 0.36)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
