\name{cjs_gamma}
\alias{cjs_gamma}
\alias{ms_gamma}
\title{HMM Transition matrix functions}
\usage{
  cjs_gamma(id,ddl,parlist,parameters)
    ms_gamma(id,ddl,parlist,parameters)
}
\arguments{
  \item{id}{sequential id for the observed sequence}

  \item{ddl}{design data list of parameters for each id}

  \item{parlist}{list of parameter vectors split by type
  (eg Phi, p in CJS)}

  \item{parameters}{formulas for each parameter type}
}
\value{
  list of occasion-specific transition matrices - Gamma in
  Zucchini and MacDonald (2009)
}
\description{
  Functions that compute the transition matrix for various
  models. Currently only CJS and MS models are included.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}
\references{
  Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov
  Models for Time Series: An Introduction using R. Chapman
  and Hall, Boca Raton, FL. 275p.
}

