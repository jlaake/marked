\name{fix.parameters}
\alias{fix.parameters}
\title{Fixing real parameters in crm models}
\usage{
  fix.parameters(x, subset = NULL, occasions, values)
}
\arguments{
  \item{x}{processed data list}

  \item{subset}{logical expression for set of animals; used
  in subset function}

  \item{occasions}{vector of occasion numbers}

  \item{values}{either single or vector of values; if
  latter, must match length of occasions}
}
\value{
  Matrix with 3 columns: 1) animal row number, 2) occasion,
  3) value for real parameter
}
\description{
  Creates matrix with appropriate structure for fixing real
  parameters in a crm model for a specified set of capture
  histories and occasions.
}
\author{
  Jeff Laake
}

