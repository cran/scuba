\name{deco.ceiling}
\alias{deco.ceiling}
\title{
  Decompression Ceiling for a Diver
}
\description{
  Given the state of a diver's tissues,
  find the diver's decompression ceiling depth (minimum permissible depth)
  or decompression ceiling pressure (minimum permissible ambient
  pressure).
}
\usage{
deco.ceiling(state, model, what = c("depth", "pressure"))
}
\arguments{
  \item{state}{
    Diver's tissue saturations. See Details.
  }
  \item{model}{
    Haldane type decompression model (object of class \code{"hm"})
    or one of the standard model names recognised by \code{\link{pickmodel}}.
  }
  \item{what}{
    What to calculate. See Details.
  }
}
\details{
  Given the current state of a diver's tissues,
  this command calculates the diver's decompression ceiling:
  the shallowest depth (or lowest ambient pressure) to which it is safe to
  ascend.

  The argument \code{state} gives the diver's current state
  or the diver's progressive state during a dive.
  It has typically been calculated using \code{\link{haldane}}.
  It should be either:
  \itemize{
    \item
    A numeric vector giving the diver's current nitrogen tensions
    for each tissue. 
    The length of the vector matches the number of tissues in
    the \code{model}.
    \item
    A matrix giving the diver's current nitrogen and (if present)
    helium tensions for each tissue. The rows of the matrix correspond to
    tissues in the model. The first column (or the column labelled
    \code{"N2"}) contains the nitrogen tensions, and if present, the second
    column (or column labelled \code{"He"}) contains the helium
    tensions.
    \item
    A three-dimensional array giving the progressive states of
    the diver over time. The first index corresponds to time.
    The second index corresponds to tissues in the model. The third
    coordinate corresponds to inert gases \code{"N2"} and \code{"He"}.
  }
  Entries in the vector, matrix or array \code{state} are pressures
  in atmospheres absolute (ata).

  The minimum tolerated ambient pressure \code{Pamb.tol} is calculated
  using the equivalent of Buehlmann's formula (Buehlmann et al, 1995,
  section 6.5, page 117). If \code{what="pressure"}, this minimum
  tolerated ambient pressure is returned (negative pressures are reset
  to zero). If \code{what="depth"} this pressure is converted to
  diving depth in metres assuming that the surface pressure is 1 ata
  (negative depths are reset to zero).
}
\value{
  A numeric vector, matrix or array of the same type as \code{state}
  and with the same dimensions as \code{state},
  giving the decompression ceiling for each tissue,
  each inert gas, and each time point as given.
  If \code{what="pressure"} the entries in the result are pressures
  in ata. If \code{what="depth"} the entries are depths in metres.
}
\references{
  Buehlmann, A.B., Voellm, E.B. and Nussberger, P. (2002)
  Tauchmedizin. 5th edition. Springer-Verlag, Berlin.
  ISBN 3-540-42979-4.
}
\seealso{
  \code{\link{hm}} for decompression models.
  
  \code{\link{haldane}} to calculate tissue tensions
  during or after a dive.
}
\examples{
   # diver state after saturation diving on EANx30 at 20 metres
   s <- saturated.state("Z", 20, nitrox(0.30))

   # decompression ceiling for each compartment
   depthceiling <- deco.ceiling(s, "Z")
   round(depthceiling, 2)
   round(max(depthceiling), 2)
   # decompress at 8 metres or deeper

   # flying-after-diving with cabin pressure 0.8 ata
   pressureceiling <- deco.ceiling(s, "Z", "pressure")
   round(pressureceiling, 2)
   any(pressureceiling > 0.8)
   # immediate flying-after-diving is not safe
}
\author{\adrian.}
\keyword{utilities}
