\name{air}
\alias{air}
\alias{is.air}
\title{Air}
\description{
  A dataset representing air as a breathing gas,
  and a function which tests whether a gas is air.
}
\usage{
  air
  is.air(g)
}
\arguments{
  \item{g}{
    An object of class \code{"gas"}.
  }
}
\details{
  An object of class \code{"gas"} represents a breathing gas.
  Such objects are required for various calculations in the
  \code{scuba} library.

  The dataset \code{air} represents compressed air (21\% oxygen,
  79\% nitrogen) as a breathing gas. It is equivalent to
  \code{nitrox(0.21)}.

  The function \code{is.air} expects its argument \code{g}
  to be a gas object. It returns \code{TRUE} if 
  \code{g} is equivalent to air, and \code{FALSE} otherwise.
}
\seealso{
  \code{\link{nitrox}}
}
\examples{
  v <- air
  is.air(v)
  v <- nitrox(0.21)
  is.air(v)
  # both are TRUE
}
\author{\adrian.}
\keyword{utilities}
