\name{make.probitCJS.data}
\alias{make.probitCJS.data}
\title{Create release-recapture history data for MCMC}
\usage{
  make.probitCJS.data(data, id.col, ch.cols,
    covariate.cols)
}
\arguments{
  \item{data}{A data frame where each row contains the
  capture history and covariates for each individual. The
  cature history must be spread over multiple columns}

  \item{id.col}{Integer or character indicating the column
  that holds a personal ID for each individual}

  \item{ch.cols}{Vector of integers giving the columns that
  contain the capture histories}

  \item{covariate.cols}{Vector of integers that indicate
  the columns containg covariate information}
}
\value{
  A data frame
}
\description{
  Creates needed constructs from the release-recapture
  history.
}
\author{
  Devin Johnson
}

