\name{probitCJS}
\alias{probitCJS}
\title{Perform MCMC analysis of a CJS model}
\usage{
  probitCJS(p.model = ~1, phi.model = ~1, data, prior.list,
  burnin, iter, init.list)
}
\arguments{
  \item{p.model}{A formula for the detection portion of the
  CJS model}

  \item{phi.model}{A formula for the survival portion of
  the CJS model}

  \item{data}{A probitCJS.data object containg the data for
  the model}

  \item{prior.list}{A named list containing \item{phi}{A
  named list containg 'mu', the prior mean and 'Q' the
  prior precison matrix for the survival parameters}
  \item{p}{A named list containg 'mu', the prior mean and
  'Q' the prior precison matrix for the survival
  parameters}}

  \item{burnin}{number of iteration to initially discard
  for MCMC burnin}

  \item{iter}{number of iteration to run the Gibbs sampler
  for following burnin}

  \item{init.list}{A named list with initial values [more
  on this later]}
}
\value{
  A list with MCMC output [I'll add to this later]
}
\description{
  Takes an object of class probitCJS.data created with the
  function \link{make.probitCJS.data} and draws a sample
  from the posterior distribution using a Gibbs sampler
}
\examples{
# Analysis of the dipper data
data(dipper)
dipper <- splitCH(data=dipper)
dipper.mcmc <- make.probitCJS.data(dipper, ch.cols=3:9, covariate.cols=2)
fit1 <- probitCJS(p.model=~time*sex, phi.model=~time*sex, data=dipper.mcmc, burnin=100, iter=1000)
# Real parameter values
fit1$reals.phi
fit1$reals.p
}
\author{
  Devin Johnson
}

