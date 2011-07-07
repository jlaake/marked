\name{print.crm}
\alias{print.crm}
\alias{coef.crm}
\title{Print model results or extract coefficients}
\description{
Provides a printed simple summary of the model results or extracts the beta
coefficients from the model results.
}
\usage{
\method{print}{crm}(x,...)
\method{coef}{crm}(object,...)
}
\arguments{
  \item{x}{crm model result}
  \item{object}{crm model result}
  \item{...}{generic arguments not used here}
}
\value{
\code{print} prints a simple summary of the model to the screen and returns NULL.
\code{coef} returns a dataframe with estimates and standard errors and confidence intervals
if hessian=TRUE on model run.}
\author{Jeff Laake }
\seealso{ \code{\link{crm}} }
\keyword{utility}


