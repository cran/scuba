\name{whichtank}
\alias{whichtank}
\alias{whichtank<-}
\title{Which Tanks are Used during a Dive}
\description{
  Determine which tank of breathing gas is used
  at each time point during a dive.
}
\usage{
whichtank(d)
whichtank(d) <- value
}
\arguments{
  \item{d}{The dive (an object of class \code{"dive"}).}
  \item{value}{
    Vector of integers or character strings identifying
    which tank of breathing gas is used at each time point.
  }
}
\details{
  An object of class \code{"dive"} represents a scuba dive,
  including information about depth, time and gas breathed at each
  stage of the dive. These objects are created by the function
  \code{\link{dive}}.

  A dive object has a \emph{tank list} which is a list of the
  tanks of breathing gas that were used (or were available to be used)
  during the dive. The function \code{\link{tanklist}} returns
  this list.

  The selection of tanks, i.e. which tank is actually used at each stage of the
  dive, is specified by \code{whichtank}. The command \code{whichtank(d)}
  returns a vector of integers or character strings,
  identifying which tank in the tank list is in use at each
  waypoint during the dive. That is, \code{whichtank(d)[i]} is the tank
  in use at the \code{i}th waypoint during the dive. 

  The command \code{whichtank(d) <- value} will change the selection of
  tanks used at each stage during the dive. Here \code{value} should be
  a vector of integers or character strings identifying tanks in the
  tank list, and the length of \code{value} should be 
  equal to the length of the vectors \code{depths.dive(d)} and
  \code{depths.dive(d)}.

  A common use of \code{whichtank} is to specify that a particular
  gas should be used only in a particular range of depths. This is done
  by applying \code{ifelse} to \code{\link{depths.dive}} as shown
  in the Examples.
}
\value{
  The value of \code{whichtank(d)} is a vector of integers or
  character strings.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{dive}},
  \code{\link{tanklist}}
}
\examples{
   # tanks are numbered
   d <- dive(air, c(30,40), 6, nitrox(0.5), c(6,3), c(3,3))
   d
   tanklist(d)
   whichtank(d)
   # change choice of tank at 6 metre deco stop
   # The Hard Way:
   whichtank(d) <- c(1,1,1,1,1,2,2,2)
   d

   # The Snazzy Way:
   # if shallower than 3 metres, then tank 2, else tank 1
   whichtank(d) <- ifelse(depths.dive(d) <= 3, 2, 1)
   d

   # tanks with names
   dd <- dive(tanklist=list(travel=trimix(0.18, 0.5), deco=nitrox(0.8)),
             tank="travel", c(30,40), c(20, 10), 9, tank="deco",
             c(9,10), c(6,5), c(3,5))
   dd
   # if shallower than 6 metres, then deco gas, else travel gas
   whichtank(dd) <- ifelse(depths.dive(dd) <= 6, "deco", "travel")
   dd
}
\keyword{utilities}
