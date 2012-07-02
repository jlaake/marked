\name{js.hessian}
\alias{js.hessian}
\title{Compute variance-covariance matrix for fitted JS model}
\usage{
  js.hessian(model, nobstot)
}
\arguments{
  \item{model}{fitted JS model from function crm}

  \item{nobstot}{number of unique caught at least once by
  group if applicable}
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

