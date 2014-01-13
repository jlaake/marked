\name{print.crm}
\alias{coef.crm}
\alias{print.crm}
\alias{print.crmlist}
\title{Print model results or extract coefficients}
\usage{
\method{coef}{crm}(object,...)
       \method{print}{crm}(x,...)
       \method{print}{crmlist}(x,...)
}
\arguments{
  \item{x}{crm model result or list of model results}

  \item{object}{crm model result}

  \item{...}{generic arguments not used here}
}
\value{
\code{print} prints a simple summary of the model to the
screen and returns NULL. \code{coef} returns a dataframe
with estimates and standard errors and confidence intervals
if hessian=TRUE on model run.
}
\description{
Provides a printed simple summary of the model results or
extracts the beta coefficients from the model results.
}
\author{
Jeff Laake
}
\seealso{
\code{\link{crm}}
}
\keyword{utility}

