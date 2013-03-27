\name{mixed.model.admb}
\alias{mixed.model}
\alias{mixed.model.admb}
\alias{mixed.model.dat}
\alias{reindex}
\title{Mixed effect model contstruction}
\usage{
  mixed.model.admb(formula,data)

  mixed.model(formula,data)

  mixed.model.dat(x,con,idonly,n)

  reindex(x,id)
}
\arguments{
  \item{formula}{formula for mixed effect mode in the form
  used in lme4; ~fixed +(re1|g1) +...+(ren|gn)}

  \item{data}{dataframe used to construct the design
  matrices from the formula}

  \item{x}{list structure created by mixed.model.admb}

  \item{con}{connection to data file which contents will be
  appended}

  \item{id}{vector of factor values used to split the data
  up by individual capture history}

  \item{idonly}{TRUE, if random effects not crossed}

  \item{n}{number of capture history records}
}
\value{
  mixed.model.admb returns a list with elements fixed.dm,
  the design matrix for the fixed effects; re.dm, a
  combined design matrix for all of the random effects; and
  re.indices, matrix of indices into a single vector of
  random effects to be applied to the design matrix
  location. mixed.model returns similar quantities for the
  random effects in a list structure (re.list) except that
  the indices are limited to the particular random effect
  grouping. May be more useful with R than ADMB.
}
\description{
  Functions that develop structures needed for a mixed
  effect model
}
\details{
  mixed.model.admb - creates design matrices and supporting
  index matrices for use of mixed model in ADMB

  mixed.model - creates design matrices and supporting
  index matrices in an alternate list format that is not as
  easily used in ADMB

  mixed.model.dat - writes to data file (con) for fixed and
  random effect stuctures

  reindex - creates indices for random effects that are
  specific to the individual capture history; it takes
  re.indices, splits them by id and creates a ragged array
  by id (used.indices) with the unique values for that id.
  index.counts is the number of indices per id to read in
  ragged array. It then changes re.indices to be an index
  to the indices within the id from 1 to the number of
  indices within the id.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}

