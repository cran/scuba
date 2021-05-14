\name{param}
\alias{param}
\title{Extract parameters from Haldane model}
\description{
  Extracts the tissue compartment parameters from a Haldane type
  decompression model.
}
\usage{
  param(model, species = "N2", what = "HalfT")
}
\arguments{
  \item{model}{
    The decompression model. An object of class \code{"hm"}.
  }
  \item{species}{
    Character string.
    The inert gas for which parameters are required:
    \code{"N2"} for nitrogen or \code{"He"} for helium.
  }
  \item{what}{
    Character string.
    The parameter required: \code{"HalfT"} for half-times,
    \code{"M0"} for surfacing M-values, or
    \code{"dM"} for slopes of M-values.
  }
}
\details{
  This command extracts the tissue compartment parameters
  from a Haldane-type decompression model.

  The argument \code{model} should be an object of class \code{"hm"}
  (created by \code{\link{hm}} or obtained using
  \code{\link{pickmodel}}).

  See \code{\link{hm}} for explanation of the parameters.
}
\value{
  A vector containing the values of the specified parameter for
  each tissue compartment. If the model does not contain this data,
  \code{NULL} is returned.
}
\seealso{
  \code{\link{hm}}, 
  \code{\link{pickmodel}}
}
\examples{
   # halftimes for the DSAT model
   param(pickmodel("DSAT"), "N2", "HalfT")
}
\author{\adrian.}
\keyword{utilities}
