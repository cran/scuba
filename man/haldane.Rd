\name{haldane}
\alias{haldane}
\title{Tissue Saturation by Haldane Model}
\description{
  Computes a diver's tissue saturation 
  during and after a dive, as predicted by a Haldane model.
}
\usage{
  haldane(d,
	 model=pickmodel("DSAT"),
	 prevstate=NULL,
         progressive=FALSE,
         relative=FALSE,
         deco=FALSE)
}
\arguments{
  \item{d}{
    The dive profile. An object of class \code{\link{dive}}. 
  }
  \item{model}{
    The decompression model. An object of class \code{"hm"}.
    Defaults to the DSAT (PADI) model.
  }
  \item{prevstate}{
    Optional. Initial state of the diver.
    A data frame containing the tissue saturations 
    for each tissue compartment in the model, at the start of the dive.
    Defaults to the state of a diver with no previous dive history.
  }
  \item{progressive}{
    Logical flag. If \code{TRUE}, the tissue saturations are computed
    at every time point during the dive. If \code{FALSE} (the default),
    only the final tissue saturation at the end of the dive is computed.
  }
  \item{relative}{
    Logical flag indicating whether to compute relative tissue
    saturations. If \code{FALSE} (the default),
    tissue saturations are expressed as absolute pressures in
    atmospheres absolute (ata). If \code{TRUE}, the tissue saturation for each
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
}
\value{
  If \code{relative=FALSE} and \code{progressive=FALSE}, a data frame giving
  the diver's inert gas saturation state at the end of the
  dive. Each row of the data frame corresponds to a tissue compartment.
  The column \code{N2} gives the nitrogen tension (in atmospheres absolute) of
  each compartment. The column \code{He}, if present, gives the
  Helium tension (in atmospheres absolute) in each compartment.

  If \code{relative=FALSE} and \code{progressive=TRUE},
  a three-dimensional array giving
  the diver's inert gas saturation state at each time point during the
  dive. The first dimension of the array corresponds to successive time points
  during the dive (the times can be extracted by
  \code{\link{times.dive}}).
  The second dimension corresponds to the tissue compartments.
  The third dimension corresponds to the inert gases
  (\code{N2} and/or \code{He}).
  The entries are gas tensions (in atmospheres absolute).

  If \code{relative=TRUE} and \code{progressive=FALSE}, a vector
  giving the diver's relative saturation of inert gas at the end of the
  dive. Entries in the vector correspond to tissue compartments.
  The entries are relative gas tensions, that is,
  the total inert gas (Nitrogen plus Helium) tissue saturation
  divided by the appropriate M-value for that compartment: either
  the surfacing M-value (if \code{deco=FALSE}) or the
  M-value at current depth (if \code{deco=TRUE}).

  If \code{relative=TRUE} and \code{progressive=TRUE},
  a matrix giving the diver's relative saturation
  at each time point during the
  dive. Rows of the array correspond to successive time points
  during the dive. Columns correspond to the tissue compartments.
  The entries are relative gas tensions, that is,
  the total inert gas (Nitrogen plus Helium) tissue saturation
  divided by the appropriate M-value for that compartment: either
  the surfacing M-value (if \code{deco=FALSE}) or the
  M-value at current depth (if \code{deco=TRUE}).
}
\details{
  This command computes a diver's nitrogen saturation 
  during and after a dive, as predicted by a Haldane model
  (Boycott et al, 1908).

  A Haldane-type decompression model describes the diver's body
  as a set of independent compartments
  connected directly to the breathing gas and governed by classical
  diffusion.

  Henry's Law is applied 
  to predict the on- and off-gassing of inert gas
  in each tissue (compartment) of the model.
  The resulting differential equations are solved analytically
  (the solution is often called the \sQuote{Schreiner equation} in
  the decompression literature).
  
  The argument \code{prevstate} represents the tissue saturation of the
  diver at the start of the dive. It should be a data frame, with one row for
  each compartment of the decompression model, and one column for each
  inert gas (\code{N2} and/or \code{He}) in the model. Such data frames
  are usually generated by \code{\link{saturated.state}} or
  \code{haldane}.

  If \code{progressive=FALSE}, the diver's tissue saturation
  at the end of the dive is calculated.

  If \code{progressive=TRUE},
  the tissue saturations are calculated at each 
  waypoint during the dive. The corresponding times
  (extracted by \code{\link{times.dive}(d)})
  are \emph{not equally spaced} over time.

  To compute the tissue saturation
  at an arbitrary instant of time during the dive, \code{tim}, use
  \code{haldane(chop.dive(d, 0, tim))}. To view the tissue saturation
  at arbitrary instants of time using interactive graphics,
  use \code{\link{showstates}}.
}
\section{Warnings}{
  Not applicable to altitude dives. Not suitable for dive planning.

  No constraints of any kind are checked. In particular
  it is not guaranteed that the model accepts the dive profile
  as a no-decompression dive.
}
\seealso{
  \code{\link{dive}},
  \code{\link{oxtox}},
  \code{\link{pickmodel}},
  \code{\link{showstates}},
  \code{\link{Mmix}},
  \code{\link{M0mix}}.
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
  
  Boycott, A.E. Damant, G.C.C. and Haldane, J.B. (1908)
  The prevention of compressed air illness.
  \emph{Journal of Hygiene} (London) \bold{8}, 342--443.

  Brubakk, A.O. and Neuman, T.S. (eds.) (2003)
  Bennett and Elliott's Physiology and Medicine of Diving.
  5th Edition. Saunders. ISBN 0-7020-2571-2

  Buehlmann, A.A. (1983)
  \emph{Dekompression - Dekompressionskrankheit}.
  Springer-Verlag.
  
  Buehlmann, A.A., Voellm, E.B. and Nussberger, P. (2002)
  \emph{Tauchmedizin}. 5e Auflage. Springer-Verlag.

  Tikvisis, P. and Gerth, W.A. (2003)
  Decompression Theory. 
  In Brubakk and Neuman (2003), Chapter 10.1, pages 419-454.

  Wienke, B.R. (1994)
  \emph{Basic diving physics and applications}.
  Best Publishing Co.
  
  Workman, R.D. (1965) Calculation of decompression schedules
  for nitrogen-oxygen and helium-oxygen dives.
  Research Report 6-65. US Navy Experimental Diving Unit.
  Washington DC. 
}
\examples{
  # First dive to 25 m for 20 min with safety stop
  d1 <- dive(c(25,20),c(5,5))
  # Evaluate saturation according to DSAT model
  s1 <- haldane(d1)
  s1
  # Look at saturation (in ata)
  barplot(s1$N2, ylab="Saturation (ata)")
  # Look at relative saturation
  M0 <- param(pickmodel("D"), "N2", "M0")
  barplot(100 * s1$N2/M0, ylab="Saturation (percent)")

  # Evaluate saturation after 2 hour surface interval
  s2 <- haldane(dive(c(0,120)), prevstate=s1)
  # Then after another dive to 18 m for 30 min with safety stop
  s3 <- haldane(dive(c(18, 30),c(5,3)), prevstate=s2)
  # Assess effect of breathing 80\% oxygen at safety stop
  s3o <- haldane(dive(c(18, 30),5, nitrox(0.8), c(5,3)), prevstate=s2)

  # Inspect saturation during dive d1 at time 10 minutes
  s10 <- haldane(chop.dive(d1, 0, 10))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}
