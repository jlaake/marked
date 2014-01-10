\name{compute.real}
\alias{compute.real}
\title{Compute estimates of real parameters}
\usage{
  compute.real(model, parameter, ddl = NULL, unique = TRUE,
    vcv = FALSE, se = FALSE, chat = 1, subset, select,
    showDesign = FALSE)
}
\arguments{
  \item{model}{model object}

  \item{parameter}{name of real parameter to be computed
  (eg "Phi" or "p")}

  \item{ddl}{list of design data}

  \item{unique}{TRUE if only unique values should be
  returned}

  \item{vcv}{logical; if TRUE, computes and returns v-c
  matrix of real estimates}

  \item{se}{logical; if TRUE, computes std errors and conf
  itervals of real estimates}

  \item{chat}{over-dispersion value}

  \item{subset}{logical expression using fields in real
  dataframe}

  \item{select}{character vector of field names in real
  that you want to include}

  \item{showDesign}{if TRUE, show design matrix instead of
  data}
}
\value{
  A data frame (\code{real}) is returned if
  \code{vcv=FALSE}; otherwise, a list is returned also
  containing vcv.real: \item{real}{ data frame containing
  estimates, and if vcv=TRUE it also contains standard
  errors and confidence intervals}
  \item{vcv.real}{variance-covariance matrix of real
  estimates}
}
\description{
  Computes real estimates and their var-cov for a
  particular subset of parameters.
}
\examples{
data(dipper)
dipper.proc=process.data(dipper,model="cjs",begin.time=1)
dipper.ddl=make.design.data(dipper.proc)
mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,
 model.parameters=list(Phi=list(formula=~sex+time),p=list(formula=~1)),hessian=TRUE)
xx=compute.real(mod.Phisex.pdot,"Phi",unique=TRUE,vcv=TRUE)
}
\author{
  Jeff Laake
}
\keyword{utility}

