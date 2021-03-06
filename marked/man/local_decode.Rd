% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMMutilities.r
\name{local_decode}
\alias{local_decode}
\title{Local decoding of HMM}
\usage{
local_decode(object, ddl = NULL, state.names = NULL)
}
\arguments{
\item{object}{fitted crm model (must be an HMM model)}

\item{ddl}{design data list}

\item{state.names}{names for states used to label output; if NULL uses strata.labels + Dead state}
}
\value{
matrix of state predictions
}
\description{
Computes state predictions one at a time for each occasion for each individual
}
\examples{
#
\donttest{
# This example is excluded from testing to reduce package check time
# cormack-jolly-seber model
data(dipper)
mod=crm(dipper,model="hmmcjs")
local_decode(mod)
}
}
\references{
Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 80.
}
\author{
Jeff Laake
}
\keyword{utility}
