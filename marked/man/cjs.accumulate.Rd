\name{cjs.accumulate}
\alias{cjs.accumulate}
\title{Accumulates common capture history values}
\usage{
  cjs.accumulate(x, model_data, nocc, freq, chunk_size)
}
\arguments{
  \item{x}{data}

  \item{model_data}{list of design matrices, fixed
  parameters and time intervals all which can vary by
  animal}

  \item{nocc}{number of capture occasions}

  \item{freq}{frequency of each capture history before
  accumulation}

  \item{chunk_size}{size that determines number of pieces
  of data/design matrix that are handled. Smaller
  chunk_size requires more time but less memory. 1e7 is
  default set in cjs.}
}
\value{
  modified model_data list that is accumulated
}
\description{
  To speed up compuation, animals with the same capture
  history and design matrix are accumulated and represented
  by a frequency. Computes starting values for Phi and p
  parameters from the list of design matrices and the
  summarized data list including ch matrix and first and
  last vectors. If any values are missing (NA) or
  abs(par)>5, they are set to 0.
}
\author{
  Jeff Laake
}
\keyword{utility}

