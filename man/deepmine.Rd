\name{deepmine}
\alias{deepmine}
\docType{data}
\title{
  Extremely Deep Decompression Dive
}
\description{
  Contains the dive profile for a long,
  extremely deep, decompression dive, using mixed gases,
  in a flooded mine.
}
\details{
  This dataset contains the dive profile data for an
  extreme dive in a flooded mine.

  This was a dive to 110 metres for 7 minutes, followed by
  147 minutes of staged decompression, including 40 minutes in
  a dry habitat at 4 metres. 
  Five different mixed gases (trimix and nitrox) were breathed.

  The diver suffered decompression sickness after this dive.
  An analysis of the dive profile was presented by Buzzacott et al (2013).

  The dataset is an object of class \code{"dive"}
  containing the dive profile.
  It can be plotted and printed (using \code{\link{plot.dive}}
  and \code{\link{print.dive}}). The nitrogen saturation can be
  computed using \code{\link{haldane}} and the
  cumulative oxygen toxicity using \code{\link{oxtox}}.
}
\format{
  An object of class \code{"dive"}.
}
\examples{
  data(deepmine)
  plot(deepmine, legendpos="right")
}
\references{
  Buzzacott, P., Papadopoulou, V., Baddeley, A., Petri, N. and Lind,
  F. (2013)
  Exceptional diving exposure with deep stops and DCS (case study).
  Conference Poster, European Underwater and Baromedical Society/
  South Pacific Underwater Medicine Society/ Southern African
  Underwater and Hyperbaric Medicine Association,
  Tricontinental Scientific Meeting on Diving
  and Hyperbaric Medicine, La Reunion, September 23-28, 2013.
}
\source{
  Peter Buzzacott.
}
\keyword{datasets}
