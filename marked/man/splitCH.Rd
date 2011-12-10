\name{splitCH}
\alias{splitCH}
\title{Split a character string of capture histories into a matrix}
\usage{
  splitCH(x = "ch", data = NULL)
}
\arguments{
  \item{x}{A vector containing the character strings of
  capture histories or the column number or name in the
  data set \code{data}}

  \item{data}{A data fram containing x if x indicates a
  column in a data frame}
}
\value{
  A data frame
}
\description{
  Split a character string of capture histories into a
  matrix. The matrix is appended to the original data set
  after
}
\author{
  Devin Johnson
}

