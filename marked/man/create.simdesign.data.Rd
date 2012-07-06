\name{create.simdesign.data}
\alias{create.simdesign.data}
\title{Create design data for simulation}
\usage{
  create.simdesign.data(num.cohorts, cohort.sizes,
    initial.age = NULL, num.occ = NULL)
}
\arguments{
  \item{num.cohorts}{number of cohorts}

  \item{cohort.sizes}{size of each cohort; if only a single
  value each cohort has same size}

  \item{initial.age}{if specified an age field is included
  as defined in details}

  \item{num.occ}{number of occasions which must be as large
  as num.cohorts; most useful for num.cohorts=1; num.occ
  >1}
}
\value{
  dataframe with columns id and design data fields: cohort
  and time and optionally age
}
\description{
  Creates a dataframe that can be used with formulae to
  specify models for parameters in a c-r simulation
  experiment.
}
\details{
  A c-r data set is composed of individuals which can be
  either caught or not caught on a series of occasions
  through time.  Each individual is part of a cohort
  defined by the time of initial capture (or release).
  Animals are assumed to age by 1 unit in the interval
  between each capture occasion from their initial age at
  first capture/release. If the length of initial.age =1,
  all animals assumed to be caught at same initial age (eg
  pups).  If it matches length of cohorts, each cohort can
  have a different initial age and if it matches number of
  animals each animal can have a different initial age.
  There are 1 or more records in the dataframe for each
  animal. For animals in the first cohort there are
  num.cohort records and for animals in the second cohort
  there are num.cohorts-1 records, etc.  As an example,
  consider 3 cohorts with 1 animal in each cohort all
  starting at age 0.The dataframe would be:

  id cohort time age 1 1 1 0 1 1 2 1 1 1 3 2 1 2 2 0 1 2 3
  1 1 3 3 0

  Added covariates can be added by including another column
  in the dataframe.  The added covariates can be specific
  to a cohort, time, age or individual.
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}

