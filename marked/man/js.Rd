\name{js}
\alias{js}
\title{Fitting function for Jolly-Seber model using Schwarz-Arnason POPAN
formulation}
\usage{
  js(x, ddl, dml, model_data = NULL, parameters,
    accumulate = TRUE, Phi = NULL, p = NULL,
    initial = NULL, method = "BFGS", hessian = FALSE,
    debug = FALSE, chunk_size = 1e+07, refit,
    itnmax = NULL, control = NULL, scale, ...)
}
\arguments{
  \item{x}{processed dataframe created by process.data}

  \item{ddl}{list of dataframes for design data; created by
  call to \code{\link{make.design.data}}}

  \item{dml}{list of design matrices created by
  \code{\link{create.dm}} from formula and design data}

  \item{model_data}{a list of all the relevant data for
  fitting the model including imat,
  Phi.dm,p.dm,Phi.fixed,p.fixed, and time.intervals. It is
  used to save values and avoid accumulation again if the
  model was re-rerun with an additional call to js when
  using autoscale or re-starting with initial values.  It
  is stored with returned model object.}

  \item{parameters}{equivalent to \code{model.parameters}
  in \code{\link{crm}}}

  \item{accumulate}{if TRUE will accumulate capture
  histories with common value and with a common design
  matrix for Phi and p to speed up execution}

  \item{Phi}{initial value of Phi; used to set intercept
  parameter}

  \item{p}{initial value of p; used to set intercept
  parameter}

  \item{initial}{initial values for parameters if desired;
  if named vector from previous run it will match to
  columns with same name}

  \item{method}{method to use for optimization; see
  \code{optimx}}

  \item{hessian}{if TRUE will compute and return the
  hessian}

  \item{debug}{if TRUE will print out information for each
  iteration}

  \item{chunk_size}{specifies amount of memory to use in
  accumulating capture histories; amount used is
  8*chunk_size/1e6 MB (default 80MB)}

  \item{refit}{non-zero entry to refit}

  \item{itnmax}{maximum number of iterations}

  \item{control}{control string for optimization functions}

  \item{scale}{vector of scale values for parameters}

  \item{...}{any remaining arguments are passed to
  additional parameters passed to \code{optimx} or
  \code{\link{js.lnl}}}
}
\value{
  The resulting value of the function is a list with the
  class of crm,js such that the generic functions print and
  coef can be used. \item{beta}{named vector of parameter
  estimates} \item{lnl}{-2*log likelihood} \item{AIC}{lnl +
  2* number of parameters} \item{convergence}{result from
  \code{optimx}; if 0 \code{optimx} thinks it converged}
  \item{count}{\code{optimx} results of number of function
  evaluations} \item{reals}{dataframe of data and real Phi
  and p estimates for each animal-occasion excluding those
  that occurred before release} \item{vcv}{var-cov matrix
  of betas if hessian=TRUE was set}
}
\description{
  A function for computing MLEs for a specified Jolly-Seber
  open population capture-recapture model for processed
  dataframe \code{x} with user specified formulas in
  \code{parameters} that create list of design matrices
  \code{dml}. This function can be called directly but is
  most easily called from \code{\link{crm}} that sets up
  needed arguments.
}
\details{
  It is easiest to call \code{js} through the function
  \code{\link{crm}}. Details are explained there.

  Be cautious with this function at present.  It does not
  include many checks to make sure values like fixed values
  will remain in the specified range of the data.  Normally
  this would not be a big problem but because
  \code{\link{js.lnl}} calls an external FORTRAN subroutine
  via \code{\link{cjs.lnl}}, if it gets a subscirpt out of
  bounds, it will cause R to terminate.  So make sure to
  save your workspace frequently if you use this function
  in its current implementation.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}
\references{
  Schwarz, C. J., and A. N. Arnason. 1996. A general
  methodology for the analysis of capture-recapture
  experiments in open populations. Biometrics 52:860-873.
}

