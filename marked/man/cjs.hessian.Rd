\name{cjs.hessian}
\alias{cjs.hessian}
\title{Compute variance-covariance matrix for fitted CJS model}
\usage{
  cjs.hessian(model, Phi.links = NULL, p.links = NULL)
}
\arguments{
  \item{model}{fitted CJS model from function crm}

  \item{Phi.links}{vector of link function names for Phi
  parameters (currently unused)}

  \item{p.links}{vector of link function names for p
  parameters (currently unused)}
}
\value{
  variance-covariance matrix for specified model
}
\description{
  A wrapper function that sets up call to hessian function
  to compute and then invert the hessian.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}

