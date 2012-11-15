\name{create.dm}
\alias{create.dm}
\alias{create.dml}
\title{Creates a design matrix for a parameter}
\usage{
  create.dm(x, formula, time.bins=NULL, cohort.bins=NULL,
    age.bins=NULL, chunk_size=1e7, remove.intercept=NULL)

  create.dml(ddl,model.parameters,design.parameters,restrict=FALSE,chunk_size=1e7)
}
\arguments{
  \item{x}{design dataframe created by
  \code{\link{create.dmdf}}}

  \item{formula}{formula for model in R format}

  \item{time.bins}{any bins of time to collapse values}

  \item{cohort.bins}{any bins of cohort to collapse values}

  \item{age.bins}{any bins of cohort to collapse values}

  \item{chunk_size}{specifies amount of memory to use in
  creating design matrices; amount used is 8*chunk_size/1e6
  MB (default 80MB)}

  \item{remove.intercept}{if TRUE, forces removal of
  intercept in design matrix}

  \item{ddl}{Design data list which contains a list element
  for each parameter type; if NULL it is created}

  \item{design.parameters}{Specification of any grouping
  variables for design data for each parameter}

  \item{model.parameters}{List of model parameter
  specifications}

  \item{restrict}{if TRUE, only use design data with Time
  >= Cohort}
}
\value{
  A design matrix constructed with the design dataframe and
  the formula.  It contains a row for each animal-occasion
  and a column for each beta parameter in the model. It
  excludes any columns that are all 0.
}
\description{
  Creates a design matrix using the design dataframe, a
  formula and any intervals defined for time, cohort and
  age.
}
\author{
  Jeff Laake
}

