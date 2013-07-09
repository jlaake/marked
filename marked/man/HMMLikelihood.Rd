\name{HMMLikelihood}
\alias{HMMLikelihood}
\alias{loglikelihood}
\title{Hidden Markov Model likelihood functions}
\usage{
  HMMLikelihood(id,x,start,m,T,fct_dmat,fct_gamma,ddl,parlist,parameters)
    loglikelihood(par,type,x,start,m,T,freq=1,fct_dmat,fct_gamma,ddl,parameters,debug=FALSE)
}
\arguments{
  \item{id}{sequential id for the observed sequence}

  \item{x}{single observed sequence (capture history) in
  HMMLikelihood and matrix of observed sequences (row:id;
  column:occasion/time) in loglikelihood}

  \item{start}{vector of 2 elements: 1) first state, 2)
  first occasion or matrix for all id's}

  \item{m}{number of states}

  \item{T}{number of occasions; sequence length}

  \item{parlist}{list of parameter vectors split by type
  (eg Phi, p in CJS)}

  \item{fct_dmat}{function to create D from parameters}

  \item{fct_gamma}{function to create gamma - transition
  matrix}

  \item{ddl}{design data list of parameters for each id}

  \item{parameters}{formulas for each parameter type}

  \item{debug}{if TRUE, print out par values and
  -log-likelihood}

  \item{par}{vector of parameter values for log-likelihood
  evaluation}

  \item{type}{vector of parameter names used to split par
  vector into types}

  \item{freq}{vector of history frequencies or 1}
}
\value{
  HMMLikelihood returns log-likelihood for a single
  sequence and loglikelihood returns the negative
  log-likelihood for all of the data.
}
\description{
  Functions HMMLikelihood and loglikelihood which compute
  the log-likelihood for a single capture-recapture
  sequence and for a set of sequences respectively. The
  function loglikelihood is called from optimizer and it in
  turn calls HMMLikelihood for each sequence in the x
  matrix. These functions are currently limited to a
  likelihood that is conditional on first capture/release
  in a known state.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}
\references{
  Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov
  Models for Time Series: An Introduction using R. Chapman
  and Hall, Boca Raton, FL. 275p.
}

