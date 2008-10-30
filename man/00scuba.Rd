\name{scuba}
\alias{scuba-package}
\alias{scuba}
\docType{package}
\title{The Scuba Package}
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
    \item gas diffusion models 
    \item decompression calculations
    \item gas toxicity calculations.
  }
  The \code{scuba} package is intended only
  for use in research and education about the mathematical
  and statistical basis of decompression theory. It is emphatically
  not designed for actual use in scuba diving and related activities.
  See the detailed disclaimer in \code{\link{scuba.disclaimer}}.
}
\section{Dive profiles}{
  A \emph{dive profile} gives the diver's depth as a function of elapsed
  time during a scuba dive. The command \code{\link{dive}} creates an
  object representing a dive profile.

  A dive profile is piecewise linear: it is a series of \emph{stages}
  that join successive \emph{waypoints}. Each waypoint is specified
  by the depth and elapsed time when it is reached. The stage between
  two waypoints is either a sojourn at a fixed depth, or an ascent or
  descent at a constant rate.

  The function \code{\link{dive}} interprets its arguments
  as a sequence of actions or events occurring during the dive.
  If an argument is a vector of length 2, it is interpreted
  as \code{c(depth,time)} specifying the depth and duration
  of a stage of the dive. If the argument is a single number,
  it is interpreted as a depth, meaning that the diver ascends or
  descends to this depth.

  For example, the command
  \code{d <- dive(c(18, 45))} specifies a dive to 18 metres
  for 45 minutes. The command \code{d <- dive(c(18, 45), c(5,3))}
  specifies a dive to 18 metres for 45 minutes followed by a safety stop
  at 5 metres for 3 minutes. Multilevel dives with any number of stages
  are specified in the same way. A dive object may include
  periods spent at the surface (depth zero) and may therefore represent
  a succession of dives separated by surface intervals. For example,
  \code{d <- dive(c(30,15),c(9,1),c(5,5),c(0,60),c(12,60),c(5,5))}
  represents two dives (with safety stops)
  separated by a one-hour surface interval.
 
  The resulting object \code{d} is an object of class \code{"dive"}.
  It can be plotted as a conventional dive profile graph by executing
  the command \code{plot(d)}. It can be printed as a table of depths and times
  by typing its name \code{d}, or by executing
  \code{print(d,seconds=FALSE)} to print times to the nearest minute.
  A summary of the dive (with such information as the average depth,
  maximum depth and the main stages of the dive) can be printed
  by typing \code{summary(d)}.

  By default, the function \code{\link{dive}} fills in some details
  about the dive. It assumes that the diver breathes compressed air;
  the dive starts and ends at the surface
  (depth zero); the diver descends 
  at the default descent rate of 30 metres per minute; and the diver
  ascends at the default ascent rate of 18 metres per minute.
  These defaults can be changed by giving extra arguments to the
  function \code{\link{dive}}.

  Dive profiles can also be modified after they are created: see below.
}
\section{Real Dive Profiles}{
  Dive profiles may also be uploaded from your dive computer
  and studied in the \pkg{scuba} package. First convert the
  uploaded profile data to a \code{data.frame} with two columns,
  the first column containing the elapsed time and the second column
  containing the depth (in metres) recorded at each time. The elapsed
  times can be either a vector of character strings
  in minutes-and-seconds format
  \code{mm:ss} or hours-minutes seconds \code{hh:mm:ss},
  or a vector of integer times measured in \emph{seconds}
  of elapsed time, or an object of class \code{\link{difftime}}
  containing the elapsed times in any time unit. Then pass this data
  frame as an argument to the function \code{\link{dive}}.

  An example of such a data frame, uploaded from a dive computer,
  is provided in the \code{\link{baron}} dataset supplied with the
  package. See the help file for
  \code{\link{baron}} for an example of how to convert a data frame
  to a \code{"dive"} object.

  The package also provides 11 real dive profiles that have already been
  converted to \code{"dive"} objects.
  See the help file for \code{\link{pedro}}.
}
\section{Decompression Calculations}{
  The \pkg{scuba} package performs the mathematical calculations of
  decompression theory:
  \itemize{
    \item the theoretical No Decompression Limit (maximum duration of a
    no-decompression dive to a specified depth) can be computed by
    \code{\link{ndl}(depth)}
    \item the quantity of nitrogen dissolved in the diver's body
    after a dive \code{d} can be computed by \code{\link{haldane}(d)}
    \item the quantity of nitrogen dissolved in the diver's body
    at each instant \bold{during} a dive \code{d} can be computed by
    \code{\link{haldane}(d, progressive=TRUE)}
    or plotted interactively by \code{\link{showstates}(d)}.
  }
  These calculations are based on the classical theory of
  decompression originated by Haldane (Boycott et al, 1908).
  The diver's body is idealised
  as a set of independent compartments, each connected directly to the
  breathing gas, and governed by classical (exponential) diffusion.

  The model parameters (the number of compartments, their diffusion
  rates, and the maximum tolerated nitrogen tension in each compartment)
  may be chosen by the user.
  By default, the model parameters are taken from the DSAT model
  which is the basis of the PADI Recreational Dive Planner.
  Alternatively, the user can choose from a variety of standard compartment
  models using the command \code{\link{pickmodel}}, or
  construct a new model using \code{\link{hm}}.
  
  No-decompression limits (the maximum duration of a no-decompression
  dive to a given depth) can be calculated using the function
  \code{\link{ndl}}. For example \code{ndl(30)} gives the 
  theoretical NDL for a dive to 30 metres,
  predicted by the DSAT model. To use the classical US Navy model instead,
  type \code{ndl(30, model="USN")} or \code{ndl(30, model=pickmodel("USN"))}.

  The nitrogen tension (the quantity of dissolved nitrogen,
  in atmospheres absolute) in the diver's body after 
  a dive, can be calculated using the function
  \code{\link{haldane}}. If \code{d} is a dive object
  then \code{haldane(d)} returns a data frame containing the
  nitrogen tissue tensions (ata) at the end of the dive,
  in each of the 8 tissue compartments
  of the DSAT model. 
  To use the US Navy model instead,
  type \code{haldane(d, "USN")} or \code{haldane(d, pickmodel("USN"))}.

  To compute the nitrogen tissue tensions at each waypoint during the
  dive, use \code{\link{haldane}(d, progressive=TRUE)}.
  
  To visualise the nitrogen tissue tensions during the dive,
  use the interactive function \code{\link{showstates}}.
  This plots the dive and waits for you to click on
  a position in the graph. The tissue tensions at that
  instant are displayed as a bar plot.

  The total oxygen toxicity incurred during a dive can be computed
  by \code{\link{oxtox}}.

  Bubble theory calculations are not yet implemented.
}
\section{Nitrox and trimix}{
  A \bold{breathing gas} is represented by an object of class \code{"gas"}.
  The object \code{air} is a representation of compressed air
  (21\% oxygen, 79\% nitrogen) as an object of this class.
  (Don't reassign another value to this object!!!)

  Nitrox mixtures (mixtures of oxygen and nitrogen) can be
  represented using the function \code{\link{nitrox}}.
  For example, EAN 32 is represented by \code{nitrox(0.32)}.

  Trimix (a mixture of oxygen, nitrogen and helium) can also be
  represented, using the command \code{\link{trimix}}.
  For example, Trimix 15/50
  (containing 15\% oxygen, 50\% helium and 35\% nitrogen)
  is represented by \code{trimix(0.15, 0.5)}.

  There are methods for \code{print} and \code{summary}
  for gas objects.

  Decompression calculations (\code{\link{haldane}},
  \code{\link{ndl}}, \code{\link{showstates}})
  also work with nitrox and trimix.

  Decompression calculations with trimix require a Haldane model
  that includes parameters for Helium diffusion.
  Use \code{pickmodel("Z")} to select the Buehlmann ZH-L16A model,
  or \code{\link{hm}} to create a new model that includes
  Helium diffusion.
  
  Standard nitrox calculations are also available, for example
  \tabular{ll}{
    \code{\link{ead}} \tab
    equivalent air depth\cr
    \code{\link{mod}} \tab
    maximum operating depth\cr
    \code{\link{maxmix}} \tab
    richest nitrox mix for a given depth
  }

  The total oxygen toxicity incurred during a nitrox or trimix dive
  can be computed by \code{\link{oxtox}}.
}
\section{Diving on different gases}{
  Every \code{"dive"} object contains information
  about the breathing gas or gases used in the dive.
  This information is determined when the \code{"dive"} object is
  created (by the function \code{\link{dive}}). The default
  breathing gas is air.
  
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

  Alternatively you can use the argument \code{tanklist} to
  specify a list of tanks of breathing gas (with optional names
  like \code{"travel"} and \code{"deco"}) and change between tanks at different
  stages of the dive using an argument of the form \code{tank=number}
  or \code{tank=name}. The tank list of a dive object can be extracted
  using \code{\link{tanklist}} and modified using
  \code{\link{tanklist}<-}.
}
\section{Manipulating dive profiles}{
  Dive profiles can also be manipulated after they are created.
  This allows you, for example, to modify the deepest portion
  of a dive (diving to a deeper depth or for a longer duration),
  to abort a dive prematurely, to cut-and-paste several dives together,
  or to consider the tissue saturation incurred by a particular
  segment of a dive.

  The commands \code{\link{depths.dive}} and \code{\link{times.dive}}
  extract the depths and elapsed times at each waypoint during the dive.
  The depths can be modified using \code{\link{depths.dive<-}}. For
  example \code{d <- dive(c(30,20))} creates a dive to 30 metres for 20
  minutes, starting and finishing at the surface;
  to change the depth to 35 metres, type \code{depths.dive(d)[2:3] <- 35}.
  Similarly the elapsed times can be modified using
  \code{\link{times.dive<-}}. It may be more convenient to use the
  functions \code{\link{durations.dive}} and
  \code{\link{durations.dive<-}} which give the duration of each stage
  (the time between two successive waypoints). 
  For example \code{durations.dive(d)[2] <- 25} would mean that the
  diver now spends 25 minutes at the bottom instead of 20 minutes.
  
  To extract only part of a dive profile, use \code{\link{chop.dive}}.

  To paste together two dive profiles or fragments of dive profiles,
  simply give them as arguments to \code{\link{dive}}.
}
\section{Changing the breathing gases}{
  A dive object has a \emph{tank list} which is a list of the
  tanks of breathing gas that were used (or were available to be used)
  during the dive. The function \code{\link{tanklist}} returns
  this list, and the function \code{\link{tanklist}<-} changes the list.

  For example,
  \code{d <- dive(c(30,20), c(5,5))} is a dive conducted using air.
  To modify it to a dive that used nitrox EANx 32, simply type
  \code{tanklist(d) <- list(nitrox(0.32))}.
  Again \code{d <- dive(air, c(30,40), 6, nitrox(0.5), c(6,3), c(3,3))} is a
  dive conducted using air (tank 1) for the deep section and EANx 50
  (tank 2) for the decompression stops at 6 metres and 3 metres.
  To change the contents of tank 1 to EANx 32, type
  \code{tanklist(d) <- list(nitrox(0.32), nitrox(0.5))} or
  just \code{tanklist(d)[[1]] <- nitrox(0.32)}. 
  To associate a name which
  each tank, give names to the entries in the tank list, for example
  \code{tanklist(d) <- list(deep=nitrox(0.32), deco=nitrox(0.5))}
  or just assign \code{names(tanklist(d)) <- c("deep", "deco")}.

  The \emph{selection} of tanks,
  i.e. which tank is actually used at each stage of the
  dive, is specified by the function \code{\link{whichtank}}.
  The command \code{whichtank(d)}
  returns a vector of integers or character strings,
  identifying which tank in the tank list is in use at each
  waypoint during the dive. That is, \code{whichtank(d)[i]} is the tank
  in use at the \code{i}th waypoint during the dive.
  The vector \code{whichtank(d)} has the same length as the vectors
  \code{depths.dive(d)} and \code{times.dive(d)}.

  To change the selection of tanks at each stage during the dive,
  use the function \code{\link{whichtank}<-}. For example,
  \code{d <- dive(air, c(30,40), 6, nitrox(0.5), c(6,3), c(3,3))} is a
  dive conducted using air (tank 1) for the deep section and EANx 50
  (tank 2) for the
  decompression stops at 6 metres and 3 metres. To change this
  so that the deco gas is only used at the 3-metre stop, type
  \code{whichtank(d) <- ifelse(depths.dive(d) < 3, 1, 2)}.
  Alternatively \code{whichtank(d)[depths.dive(d) > 3] <- 1}
  would select tank 1 for all parts of the dive deeper than 3 metres.
  These manipulations are usually easier to understand if the tanks
  have names. For example typing
  \code{names(tanklist(d)) <- c("deep", "deco")} we could then
  type \code{whichtank(d) <- ifelse(depths.dive(d) < 3, "deep", "deco")}
  or \code{whichtank(d)[depths.dive(d) > 3] <- "deep"}.
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\references{
  Bookspan, J. (1995) Diving physiology in plain English.
  Undersea and Hyperbaric Medicine Society, Kensington, Maryland (USA).
  ISBN 0-930406-13-3.
  
  Boycott, A.E. Damant, G.C.C. and Haldane, J.S. (1908)
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
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  with contributions from Vittorio Broglio
  and Pedro Antonio Neves.
}
\keyword{utilities}

