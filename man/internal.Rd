\name{scuba-internal}
\alias{M0mix}
\alias{allspecies}
\alias{as.character.gas}
\alias{as.data.frame.hm}
\alias{capable}
\alias{conform}
\alias{dive.segment}
\alias{is.dive}
\alias{is.gas}
\alias{is.rate}
\alias{gasnames}
\alias{print.gas}
\alias{print.hm}
\alias{print.rate}
\alias{print.summary.dive}
\alias{print.summary.gas}
\alias{print.summary.hm}
\alias{rate}
\alias{reconcile.df}
\alias{resolve.defaults}
\alias{summary.dive}
\alias{summary.gas}
\alias{summary.hm}
\alias{timetaken}
\title{Internal scuba functions}
\description{
  Internal scuba functions.
}
\usage{
M0mix(model, fN2, fHe)
allspecies(d, inert)
\method{as.character}{gas}(x, \dots)
\method{as.data.frame}{hm}(x, \dots)
capable(model, g, what)
conform(state, model)
dive.segment(d, i)
is.dive(x)
is.gas(x)
is.rate(x)
gasnames(fO2, fN2=1-fO2, fHe=1 - fO2 - fN2, ..., full=FALSE)
\method{print}{gas}(x, \dots)
\method{print}{hm}(x, \dots)
\method{print}{rate}(x, \dots)
\method{print}{summary.dive}(x, \dots)
\method{print}{summary.gas}(x, \dots)
\method{print}{summary.hm}(x, \dots)
rate(speed, time, up)
reconcile.df(df, tanks)
resolve.defaults(\dots)
\method{summary}{dive}(object, \dots)
\method{summary}{gas}(object, \dots)
\method{summary}{hm}(object, \dots)
timetaken(start, finish, uprate, downrate)
}

\details{
  Not yet documented.
}
\keyword{internal}
