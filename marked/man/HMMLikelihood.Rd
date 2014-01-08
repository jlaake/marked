\name{HMMLikelihood}
\alias{HMMLikelihood}
\alias{hmm.lnl}
\alias{reals}
\title{Hidden Markov Model likelihood functions}
\usage{
HMMLikelihood(par,type,x,start,m,T,freq=1,fct_dmat,fct_gamma,fct_delta,ddl,
                         dml,parameters,debug=FALSE)
       reals(ddl,dml,parameters,parlist)
       hmm.lnl(x,start,m,T,dmat,gamma,delta,freq)
}
\arguments{
  \item{x}{single observed sequence (capture history) in
  HMMLikelihood and matrix of observed sequences (row:id;
  column:occasion/time) in loglikelihood}

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

  \item{debug}{if TRUE, print out par values and
  -log-likelihood}

  \item{parlist}{list of parameter strings used to split
  par vector}

  \item{start}{for each ch, the first non-zero x value and
  the occasion of the first non-zero value}
}
\value{
HMMLikelihood returns log-likelihood for a single sequence
and loglikelihood returns the negative log-likelihood for
all of the data. reals returns either the column dimension
of design matrix for parameter or the real parameter vector
}
\description{
Function HMMLikelihood computes the log-likelihood via
hmm.like which is a wrapper for the FORTRAN code
hmm_like.f.  The function HMMlikelihood is called from
optimizer and it in turn calls HMMLikelihood for each
sequence in the x matrix.
}
\details{
Below is the deprecated R code for computation of the HMM
likelihood # R version of FORTRAN code hmmlike.f for a
single sequence
HMMLikelihood=function(x,first,m,T,dmat,gamma,delta) { #
Arguments: # x: observed sequence (capture (encounter)
history) # first: occasion to start sequence # m: number of
states # T: number of occasions; sequence length # dmat:
array of occasion specific observation probabilty matrices
# gamma: array of occasion specific transition matrices #
delta: initial state distribution # Other variables: # lnl:
log likelihood value # phi: alpha/sum(alpha) sequence as
defined in Zucchini/MacDonald # v: temp variable to hold
phi calculations # u: sum(v) # Assign prob state vector for
initial observation: delta*p(x_first)
v=delta%*%diag(dmat[first,x[first],]) # Compute
log-likelihood contribution for first observation; for #
models that condition on first observation u=1,lnl=0
u=sum(v) phi=v/u lnl=log(u) # Loop over occasions for this
encounter history (x) for(t in (first+1):T) { # Compute
likelihood contribution for this occasion
v=phi%*%gamma[t-1,,]%*%diag(dmat[t,x[t],]) u=sum(v)
lnl=lnl+log(u) # Compute updated state vector phi=v/u }
return(lnl) } R code that called R version of HMMlikelihood
for each sequence # loop over each encounter history in
sapply and # create log-likelihood vector - an element for
each x # sum is total log-likelihood across individuals #
return negative log-likelihood
neglnl=-sum(freq*sapply(1:nrow(x),function(id)
HMMLikelihood(x[id,],start[id,2],m,T,
dmat=dmat[id,,,],gamma=gamma[id,,,], delta=delta[id,])))
}
\author{
Jeff Laake <jeff.laake@noaa.gov>
}
\references{
Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models
for Time Series: An Introduction using R. Chapman and Hall,
Boca Raton, FL. 275p.
}

