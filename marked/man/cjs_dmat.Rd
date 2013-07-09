\name{cjs_dmat}
\alias{cjs_dmat}
\alias{ms_dmat}
\alias{ums_dmat}
\title{HMM Observation Probability matrix functions}
\usage{
  cjs_dmat(id,ddl,parlist,parameters)
    ms_dmat(id,ddl,parlist,parameters)
    ums_dmat(id,ddl,parlist,parameters)
}
\arguments{
  \item{id}{sequential id for the observed sequence}

  \item{ddl}{design data list of parameters for each id}

  \item{parlist}{list of parameter vectors split by type
  (eg Phi, p in CJS)}

  \item{parameters}{formulas for each parameter type}
}
\value{
  list of occasion-specific observation probability
  matrices - state-dependent distributions in Zucchini and
  MacDonald (2009)
}
\description{
  Functions that compute the probability matrix of the
  observations given the state for various models.
  Currently only CJS, MS models and MS models with state
  uncertainty are included.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}
\references{
  Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov
  Models for Time Series: An Introduction using R. Chapman
  and Hall, Boca Raton, FL. 275p.
}

