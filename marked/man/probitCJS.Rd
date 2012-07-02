\name{probitCJS}
\alias{probitCJS}
\title{Perform MCMC analysis of a CJS model}
\usage{
  probitCJS(data,
    parameters = list(Phi = list(formula = ~1), p = list(formula = ~1)),
    burnin, iter, init.list = NULL)
}
\arguments{
  \item{data}{A list of design data from make.design.data
  for model="probitCJS"}

  \item{parameters}{A model specification list with a list
  for Phi and p containing a formula and optionally a prior
  specification which is a named list containing 'mu', the
  prior mean and 'Q' the prior precision matrix}

  \item{burnin}{number of iteration to initially discard
  for MCMC burnin}

  \item{iter}{number of iteration to run the Gibbs sampler
  for following burnin}

  \item{init.list}{A named list with initial values [more
  on this later]}
}
\value{
  A list with MCMC iterations and summarized output:
  \item{beta.mcmc}{list with elements Phi and p which
  contain MCMC iterations for each beta parameter}
  \item{real.mcmc}{list with elements Phi and p which
  contain MCMC iterations for each real parameter}
  \item{beta}{list with elements Phi and p which contain
  summary of MCMC iterations for each beta parameter}
  \item{reals}{list with elements Phi and p which contain
  summary of MCMC iterations for each real parameter}
}
\description{
  Takes design data list created with the function
  \link{make.design.data} for model "probitCJS" and draws a
  sample from the posterior distribution using a Gibbs
  sampler.
}
\examples{
# Analysis of the dipper data
data(dipper)
# following example uses unrealistically low values for burnin and iteration to reduce package testing time
fit1 <- crm(dipper,model="probitCJS",model.parameters=list(Phi=list(formula=~time*sex),p=list(formula=~time+sex)), burnin=50, iter=250)
fit1
# Real parameter summary
fit1$reals
}
\author{
  Devin Johnson
}

