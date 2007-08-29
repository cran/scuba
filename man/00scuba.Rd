\name{scuba}
\alias{scuba}
\title{
The Scuba Package
}
\description{
  This is a summary of the features of 
  \code{scuba}, a package in \code{R}
  that performs theoretical calculations
  about scuba diving --- dive profiles,
  decompression models, gas toxicity and so on.
}
\details{
  \code{scuba} is a package for performing
  calculations in the theory of scuba diving.
  The package supports
  \itemize{
    \item creation, manipulation and plotting of dive profiles
    \item decompression models 
    \item gas toxicity calculations.
  }
  The \code{scuba} package is intended only
  for use in research and education about the mathematical
  and statistical basis of decompression theory. It is
  not designed for actual use in scuba diving and related activities.
  See the disclaimer in \code{\link{scuba.disclaimer}}.
}
\section{Introductory Tutorial}{
  A \bold{dive profile} (depth/time profile) is created with the command
  \code{\link{dive}}. For example, the command
  \code{d <- dive(c(18, 45))} specifies a dive to 18 metres
  for 45 minutes. The command \code{d <- dive(c(18, 45), c(5,3))}
  specifies a dive to 18 metres for 45 minutes followed by a safety stop
  at 5 metres for 3 minutes. Multilevel dives with any number of stages
  are specified in the same way.

  The resulting object \code{d} is an object of class \code{"dive"}.
  It can be plotted as a conventional dive profile by typing
  \code{plot(d)} as explained in the help entry for
  \code{\link{plot.dive}}.
  It can be printed as a table of depths and times
  by typing its name \code{d}, or by typing
  \code{print(d,seconds=FALSE)} to print times to the nearest minute,
  as explained in the help entry for \code{\link{print.dive}}.
  A summary of the dive (with such information as the average depth,
  maximum depth and the main stages of the dive) can be printed
  by typing \code{summary(d)}.

  By default, the function \code{dive} fills in some details
  about the dive. It assumes that the diver breathes compressed air;
  the dive starts and ends at the surface
  (depth zero); the diver descends 
  at the default descent rate of 30 metres per minute; and the diver
  ascends at the default ascent rate of 18 metres per minute.
  These defaults can be changed, as explained below.
}
\section{Decompression Models}{
  The quantity of nitrogen dissolved in the diver's body after a dive,
  as predicted by the classical Haldane theory, can be calculated
  using the function \code{\link{haldane}}. If \code{d} is a dive object
  then \code{haldane(d)} returns a vector containing the nitrogen
  tissue tensions (ata) at the end of the dive,
  in each of the 8 tissue compartments
  of the DSAT model. You can choose other compartment models;
  for example to use the classical US Navy model,
  type \code{haldane(d, Halftimes[["USN"]])}. For other
  alternative compartment models,
  see \code{\link{Halftimes}}, \code{\link{Workman65}}
  or \code{\link{BuehlmannL16A}}.

  To inspect the tissue tensions during the dive,
  use the interactive function \code{\link{showstates}}.
  This plots the dive and waits for you to click on
  a position in the graph. The tissue tensions at that
  instant are displayed as a bar plot.

  Bubble theory calculations are not yet implemented.
  
  The total oxygen toxicity incurred during a dive can be computed
  by \code{\link{oxtox}}.
}
\section{Nitrox calculations}{
  A \bold{breathing gas} is represented by an object of class \code{"gas"}.
  The object \code{air} is a representation of compressed air
  (21\% oxygen, 79\% nitrogen) as an object of this class.
  (Don't reassign another value to this object!!!)
  Nitrox mixtures (mixtures of oxygen and nitrogen) can be
  represented using the function \code{\link{nitrox}}.
  For example, EAN 32 is represented by \code{nitrox(0.32)}.

  There are methods for \code{print} and \code{summary}
  for gas objects.
  
  Standard nitrox calculations are also available, for example
  \tabular{ll}{
    \code{\link{ead}} \tab
    equivalent air depth\cr
    \code{\link{mod}} \tab
    maximum operating depth\cr
    \code{\link{maxmix}} \tab
    richest nitrox mix for a given depth
  }
}
\section{Diving on different gases}{
  Remember that a \code{"dive"} object contains information
  about the breathing gas used. The breathing gas is determined
  when the dive object is created, and cannot be changed later.
  Information about the breathing gas is passed to 
  the function \code{\link{dive}} when it creates the dive object.
  
  The function \code{\link{dive}} interprets its arguments
  as a sequence of actions or events occurring during the dive.
  If an argument is a vector of length 2, it is interpreted
  as \code{c(depth,time)} specifying the depth and duration
  of a stage of the dive. If the argument is a single number,
  it is interpreted as a depth, meaning that the diver ascends or
  descends to this depth.

  Each argument to \code{\link{dive}} may also be a "gas" object,
  like \code{nitrox(0.32)}, 
  which means that the diver switches to this gas.
  
  So, for example,
  \code{dive(nitrox(0.32), c(30,20))}
  means a dive to 30 metres for 20 minutes conducted on EAN 32
  (Nitrox 0.32) from start to finish.
  The command
  \code{dive(c(30,20), 5, nitrox(0.36), c(5,3))} means
  a dive on air to 30 metres for 20 minutes, ascending to 5 metres
  while breathing air, then switching to EAN 36 for a safety stop
  at 5 metres for 3 minutes.

  The total oxygen toxicity incurred by a diver on
  any Nitrox dive can be calculated using \code{\link{oxtox}}.
}
\section{Uploaded Dive Profiles}{
  Dive profiles may also be uploaded from your dive computer
  and studied in the \pkg{scuba} package. First convert the
  uploaded profile data to a \code{data.frame} with two columns,
  the first column containing the elapsed time and the second column
  containing the depth (in metres) recorded at each time. The elapsed
  times can be either a character vector in minutes-and-seconds format
  \code{mm:ss}, or a vector of integer times measured in \emph{seconds}
  of elapsed time, or an object of class \code{\link{difftime}}
  containing the elapsed times in any time unit. Then pass this data
  frame as an argument to the function \code{\link{dive}}.

  An example of a real dive profile is provided in the
  \code{\link{baron}} dataset.
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{utilities}

