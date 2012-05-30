\name{bestdoubledive}
\alias{bestdoubledive}
\title{
  Find the Best Double Dive To Given Depths
}
\description{
  For two no-decompression dives,
  to depths \code{d1} and \code{d2} metres respectively,
  separated by a given time interval, find the
  optimal durations of the two dives.
}
\usage{
bestdoubledive(d1, d2, surface = 30, verbose = FALSE, model = "D", asdive=TRUE)
}
\arguments{
  \item{d1}{
    Depth of first dive in metres.
  }
  \item{d2}{
    Depth of second dive in metres.
  }
  \item{surface}{
    Length of surface interval between dives, in minutes.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, print lots of extra information,
    and return extra information. If \code{FALSE} (the default),
    just return the optimal dive.
  }
  \item{model}{
    The decompression model. Either
    a character string, containing the name of a decompression model
    recognised by \code{\link{pickmodel}}, or 
    an object of class \code{"hm"}
    (created by \code{\link{hm}})
    representing a decompression model.
    Defaults to the DSAT model.
  }
  \item{asdive}{
    Logical. If \code{TRUE} (the default),
    the data for the optimal dive are converted into a dive object
    (object of class \code{"dive"}). If \code{FALSE}, the data are
    returned as a data frame.
  }
}
\details{
  This command implements the algorithm described by
  Baddeley and Bassom (2012) which calculates the \sQuote{best}
  double dive to two given depths
  separated by a given surface interval,
  without exceeding the no-decompression limits.

  Consider a no-decompression dive to depth \eqn{d_1}{d1} metres
  for \eqn{t_1}{t1} minutes, followed by a surface interval of \eqn{s}
  minutes, followed by a no-decompression dive to depth \eqn{d_2}{d2}
  for \eqn{t_2}{t2} minutes. 
  The \sQuote{best} double dive is defined by Baddeley and Bassom (2012)
  as the one which maximises the integral of depth
  \eqn{\Phi = t_1 d_1 + t_2 d_2}{Phi = t1 * d1 + t2 * d2}.
}
\value{
  By default (when \code{verbose=FALSE} and \code{asdive=TRUE})
  the result is a dive object (object of class \code{"dive"}).
  
  Otherwise, the result is a data frame with
  columns \code{t1} and \code{t2} containing the
  dive durations in minutes, \code{phi} containing the value of
  \eqn{\Phi}{Phi}, and \code{case} specifying which of the cases
  specified in Baddeley and Bassom (2012) provided the optimum.
  If \code{verbose=FALSE} this data frame has only one
  row, giving the best double dive. If \code{verbose=TRUE} then the
  data frame has several rows giving the candidates for optimal dive
  in each case of the algorithm. 
}
\references{
  Baddeley, A. and Bassom, A.P. (2012)
  Classical theory of decompression and the design of scuba diving tables.
  To appear in \emph{The Mathematical Scientist}.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{pickmodel}}, 
  \code{\link{hm}}, 
  \code{\link{ndl}}, 
  \code{\link{haldane}}
}
\examples{
   d <- bestdoubledive(40, 12, 15)
   plot(d)
   # Table 3 in Baddeley and Bassom (2012)
   bestdoubledive(40, 12, 15, verbose=TRUE)
   # Data for optimal dive
   bestdoubledive(40, 12, 15, asdive=FALSE)
}
\keyword{utilities}
