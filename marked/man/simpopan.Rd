\name{simpopan}
\alias{simpopan}
\title{Simulates data from POPAN model}
\usage{
  simpopan(num.cohorts = 1, Nstar = NULL,
    cohort.sizes = NULL, Phi = list(), p = list(),
    pent = list(), design.data = NULL, outfile = NULL)
}
\arguments{
  \item{num.cohorts}{number of cohorts; design is square
  with same number of c-r events as num.cohorts; number of
  recapture events is num.cohorts-1}

  \item{cohort.sizes}{a scalar giving constant size of each
  cohort or a vector of sizes of length num.cohorts}

  \item{Nstar}{a scalar giving super population size; only
  needed if cohort.sizes, design.data with cohort field, or
  pent=vector are not specified}

  \item{Phi}{a list defining the survival model with the
  following elements (see details) par - a vector of
  parameter values formula - a formula to use with
  design.data to construct model link - link function used
  with model to create probabilites (not used at present)}

  \item{p}{a list defining the capture probability model
  (same structure as Phi)}

  \item{pent}{a list defining the entry probability model
  (same structure as Phi)}

  \item{design.data}{a dataframe with design data that
  allows model construction for probabilities (see
  details).}

  \item{outfile}{prefix name of the output file for the ch
  data. extension .inp is always added for MARK}
}
\value{
  dataframe with ch (capture history) as character
}
\description{
  Simulates capture-recapture data from POPAN model with
  multiple cohorts from a super-population with defined
  entry into the population.
}
\details{
  Allows general model for Phi,p, pent or constant or
  time-specific model depending on input arguments
  Phi,p,pent. Optionally outputs file to outfile with .inp
  extension for input to program MARK. The entry of
  individuals into the population can be specified in
  several ways: 1) design.data contains cohort
  specification for each individual which is the occasion
  that it entered but not the first time it was seen. This
  allows any specification of entry into the population
  because it is done outside of the this code. 2)
  cohort.sizes are specified which are the number that
  enter for each cohort; can differ across cohorts. For
  example,
  cohort.sizes=rmultinom(1,500,prob=c(.1,.3,.2,.3,.1)) 3)
  Nstar is specifed and pent is specified as a vector; #2
  can be also done as Nstar=500,pent=c(.1,.3,.2,.3,.1) 4)
  Nstar and a formula for pent is specified; if no
  design.data are specified then default values are created
  with cohort=1; after generating an entry time for each
  animal the a cohort value is assigned to the design.data.
  If design.data are specified, all cohort values should be
  1 and the code will assign it from random entry as
  specified by formula. The pent parameters are computed
  using the mlogit transformation from the specified beta
  parameters and formula.

  The design.data dataframe contains a record for each
  animal for each time period it is contained in the study.
  The default fieldnames are id(unique # in cohort), cohort
  (cohort number from 1:num.cohorts), time (i:num.cohorts,
  where i is cohort # or all 1) and age (assumed to
  increment by 1 each occasion and defaults with 0 initial
  age) The function create.simdesign.data can be used to
  create an initial dataframe which can be supplemented
  with other covariates. If a formula is given but no
  design.data is not provided this function calls
  create.simdesign.data to construct a default dataframe
  for the specified problem.

  For designation of models for Phi or p, there are 3
  options here: constant model: Phi=list(par=value) time
  model: Phi=list(par=c(val1,val2,...valk)) where
  k=num.cohorts-1 is number of survival intervals general
  model:
  Phi=list(par=c(val1,val2,...valk),formula=~yourformula))
  k is number of cols in model matrix

  For first 2, identity link is assumed but logit link can
  be used. For formula, logit link is assumed and required.
  The formula must match fields used in design.data.  Note
  that fields in design.data are numeric and the formula
  must specify as.factor to create a factor variable.  For
  example a time model with 3 times and the same survivals
  in each can be specified as:

  Phi=list(par=c(0.731,0.786)) or
  s=list(par=c(1,0.3),formula=~as.factor(time))

  The above are not equivalent to
  Phi=list(par=c(1,0.3),formula=~time) which is where
  Phi=exp(1+0.3*time)/(1+exp(1+0.3*time)), Phi1=0.786,
  Phi2=0.831.  Note that survivals are indexed by time at
  beginning of interval 1 for interal 1 to 2 and capture
  probabilities by time of occasion.  Since these are truly
  recapture probabilities it starts with occasion 2, so
  p=list(par=c(1,0.3),formula=~time) would give values:
  p2=0.832, p3=0.870
}
\examples{
library(RMark)
do.popan <- function(Phi,p,cohort.sizes,reps,...)
#
#  do.popan -  a simple example piece of code to show simulation/analysis loop
#              with a POPAN model
#
#  Arguments:
#
#  Phi          -  parameter list for Phi
#  p            -  parameter list for p
#  cohort.sizes -  fixed entry cohort.sizes
#  reps         -  number of simulation reps
#  ...          -  arguments passed to mark
#
#
#  Value:
#
#  results - for this simple example it will be a matrix of the real parameter estimates
#
#  Functions used: simpopan, mark
#
#
{
	results=NULL
	for(i in 1:reps)
	{
		cat("rep ",i,"\\n")
		xx=simpopan(length(cohort.sizes),cohort.sizes,Phi=Phi,p=p,pent=NULL)
		mod<-mark(xx,model="POPAN",title="sim test",output=FALSE,...)
		results=rbind(results,mod$results$real$estimate)
	}
	return(results)
}
library(marked)
do.js.compare=function(Phi,p,num.cohorts,cohort.sizes,reps,...)
#
#  do.js.compare -  a simple example piece of code to show simulation/analysis loop
#                   with a JS model and analysis by RMark/MARK and marked
#
#  Arguments:
#
#  Phi          parameter list for Phi
#  p            parameter list for p
#  num.cohorts  number of cohorts
#  cohort.sizes size or sizes of cohorts
#  reps         number of simulation reps
#  ...          additional optional arguments passed to mark
#
#
#  Value:
#
#  results - for this simple example it will be a matrix of the beta parameter estimates and
#	         a matrix of the standard errors
#
#  Functions used: simcjs, mark
#
#
{
	results=NULL
	results.se=NULL
	crmresults=NULL
	crmresults.se=NULL
	for(i in 1:reps)
	{
		cat("rep ",i,"\\n")
		xx=simpopan(num.cohorts,cohort.sizes=cohort.sizes,Phi=Phi,p=p)
		mod<-mark(xx,model="POPAN",output=FALSE,...)
		results=rbind(results,mod$results$beta$estimate)
		results.se=rbind(results.se,mod$results$beta$se)
		mod<-crm(xx,model="js",hessian=TRUE)
		crmresults=rbind(crmresults,mod$beta)
		crmresults.se=rbind(crmresults.se, sqrt(diag(mod$vcv)))
		gc()
	}
	return(list(mark=list(beta=results,se=results.se),crm=list(beta=crmresults,se=crmresults.se)))
}
xx=do.js.compare(.9,.4,7,100,25)
summary(xx$mark$beta-xx$crm$beta)
summary(xx$mark$se-xx$crm$se)
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}

