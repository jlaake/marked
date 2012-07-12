\name{cjs.initial}
\alias{cjs.initial}
\title{Computes starting values for CJS p and Phi parameters}
\usage{
  cjs.initial(dml, imat, link = "logit")
}
\arguments{
  \item{dml}{design matrix list for Phi and p}

  \item{imat}{list containing chmat, first and last}

  \item{link}{either "logit" (for cjs) or "probit" (for
  probitCJS)}
}
\value{
  vector of initial parameter estimates
}
\description{
  Computes starting values for Phi and p parameters from
  the list of design matrices and the summarized data list
  including ch matrix and first and last vectors. If any
  values are missing (NA) or abs(par)>5, they are set to 0.
}
\author{
  Jeff Laake
}
\keyword{utility}

