\name{probitMsCJS}
\alias{probitMsCJS}
\title{Perform MCMC analysis of a CJS model}
\usage{
probitMsCJS(ddl, dml, parameters, design.parameters, burnin, iter,
  initial = NULL, imat = NULL)
}
\arguments{
  \item{ddl}{list of dataframes for design data; created by
  call to \code{\link{make.design.data}}}

  \item{dml}{list of design matrices created by
  \code{\link{create.dm}} from formula and design data}

  \item{parameters}{A model specification list with a list
  for Phi and p containing a formula and optionally a prior
  specification which is a named list containing 'mu', the
  prior mean and 'tau' the scale for the conjugate prior
  precision matrix X'X.}

  \item{design.parameters}{Specification of any grouping
  variables for design data for each parameter}

  \item{burnin}{number of iteration to initially discard
  for MCMC burnin}

  \item{iter}{number of iteration to run the Gibbs sampler
  for following burnin}

  \item{initial}{A named list (Phi,p). If null and imat is
  not null, uses cjs.initial to create initial values;
  otherwise assigns 0}

  \item{imat}{A list of vectors and matrices constructed by
  \code{\link{process.ch}} from the capture history data}
}
\value{
A list with MCMC iterations and summarized output:
\item{beta.mcmc}{list with elements Phi and p which contain
MCMC iterations for each beta parameter}
\item{real.mcmc}{list with elements Phi and p which contain
MCMC iterations for each real parameter} \item{beta}{list
with elements Phi and p which contain summary of MCMC
iterations for each beta parameter} \item{reals}{list with
elements Phi and p which contain summary of MCMC iterations
for each real parameter}
}
\description{
Takes design data list created with the function
\link{make.design.data} for model "probitMsCJS" and draws a
sample from the posterior distribution using a Gibbs
sampler.
}
\examples{
\donttest{
# Analysis of the dipper data
data(dipper)
# following example uses unrealistically low values for burnin and iteration
# to reduce package testing time
fit1 <- crm(dipper,model="probitCJS",model.parameters=list(Phi=list(formula=~time*sex),
p=list(formula=~time+sex)), burnin=100, iter=1000)
fit1
# Real parameter summary
fit1$results$reals
}
}
\author{
Devin Johnson
}

