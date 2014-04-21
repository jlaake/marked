\name{R_HMMLikelihood}
\alias{R_HMMLikelihood}
\alias{backward_prob}
\alias{loglikelihood}
\title{Hidden Markov Functions}
\usage{
R_HMMLikelihood(x,first,m,T,dmat,gamma,delta)
       loglikelihood(par,type,x,id,start,m,T,freq=1,fct_dmat,fct_gamma,fct_delta,
                        ddl,dml,parameters)
       backward_prob(x,first,m,T,dmat,gamma)
}
\arguments{
  \item{x}{single observed sequence (capture history)}

  \item{id}{id of capture history for computation}

  \item{first}{occasion to initiate likelihood calculation
  for sequence}

  \item{m}{number of states}

  \item{T}{number of occasions; sequence length}

  \item{dmat}{observation probability matrices}

  \item{gamma}{transition matrices}

  \item{delta}{initial distribution}

  \item{par}{vector of parameter values for log-likelihood
  evaluation}

  \item{type}{vector of parameter names used to split par
  vector into types}

  \item{freq}{vector of history frequencies or 1}

  \item{fct_dmat}{function to create D from parameters}

  \item{fct_gamma}{function to create gamma - transition
  matrix}

  \item{fct_delta}{function to create initial state
  distribution}

  \item{ddl}{design data list of parameters for each id}

  \item{dml}{list of design matrices; one entry for each
  parameter; each entry contains fe and re for fixed and
  random effects}

  \item{parameters}{formulas for each parameter type}

  \item{start}{for each ch, the first non-zero x value and
  the occasion of the first non-zero value}
}
\value{
both return log-likelihood, alpha and phi vectors and dmat
and gamma matrices for a single sequence
}
\description{
R implementation of HMMs described in processed report.
Function HMMLikelihood renamed to R_HMMLikelihood and
loglikelihood modified to work with a single x and to
return lnl, alpha, phi, v, dmat, and gamma values.
backward_prob returns the backward probabilities. These are
not used by the fitting code.
}
\author{
Jeff Laake <jeff.laake@noaa.gov>
}
\references{
Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models
for Time Series: An Introduction using R. Chapman and Hall,
Boca Raton, FL. 275p.
}

