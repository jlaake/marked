\name{simcjs}
\alias{simcjs}
\title{Simulates data from CJS model}
\usage{
  simcjs(num.cohorts = 1, cohort.sizes, Phi = list(),
    p = list(), design.data = NULL, outfile = NULL)
}
\arguments{
  \item{num.cohorts}{number of cohorts; design is square
  with same number of c-r eventsas num.cohorts; number of
  recapture events is num.cohorts-1; it does not return ch
  values for last release cohort}

  \item{cohort.sizes}{a scalar giving constant size of each
  cohort or a vector of sizes of length num.cohorts}

  \item{Phi}{a list defining the survival model with the
  following elements (see details) par - a vector of
  parameter values formula - a formula to use with
  design.data to construct model link - link function used
  with model to create probabilites (not used at present)}

  \item{p}{a list defining the capture probability model
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
  Simulates capture-recapture data from CJS model of
  releases of multiple cohorts. Allows general model for
  Phi and p or constant or time-specific model depending on
  input arguments Phi and p. Optionally outputs file to
  outfile with .inp extension for input to program MARK.
}
\details{
  The design.data dataframe contains a record for each
  animal for each time period it is contained in the study.
  The default fieldnames are id(unique # in cohort), cohort
  (cohort number from 1:num.cohorts), time (i:num.cohorts,
  where i is cohort #) and age (assumed to increment by 1
  each occasion and defaults with 0 initial age) The
  function create.simdesign.data can be used to create an
  initial dataframe which can be supplemented with other
  covariates. If a formula is given but no design.data is
  not provided this function calls create.simdesign.data to
  construct a default dataframe for the specified problem.

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
do.cjs <- function(Phi,p,reps,...)
#
#  do.cjs -  a simple example piece of code to show simulation/analysis loop
#            with a CJS model
#
#  Arguments:
#
#  Phi -  parameter list for Phi
#  p   -  parameter list for p
#  reps-  number of simulation reps
#  ... -  arguments passed to mark
#
#
#  Value:
#
#  results - for this simple example it will be a matrix of the real parameter estimates
#
#  Functions used: simcjs, mark
#
#
{
	results=NULL
	for(i in 1:reps)
	{
		cat("rep ",i,"\\n")
		xx=simcjs(3,500,Phi=Phi,p=p)
		mod<-mark(xx,title="sim test",output=FALSE,...)
		results=rbind(results,mod$results$real$estimate)
	}
	return(results)
}
do.cjscov <-	function(Phi,p,reps)
#
#  do.cjscov -  a simple example piece of code to show simulation/analysis loop
#               with a CJS model including a single covariate
#
#  Arguments:
#
#  Phi -  parameter list for Phi
#  p   -  parameter list for p
#  reps-  number of simulation reps
#
#
#  Value:
#
#  results - for this simple example it will be a list with 2 elements: a matrix of the real parameter estimates
#            and a matrix of beta estimates
#
#  Functions used: create.simdesign.data, simcjs, mark
#
#
{
	results.real=NULL
	results.beta=NULL
	for(i in 1:reps)
	{
		cat("rep ",i,"\\n")
		cov=data.frame(id=1:1500,cov=rnorm(1500,0,1))
		ddi=create.simdesign.data(3,500,0)
		ddi=merge(cov,ddi,by="id",sort=FALSE)
		xx=simcjs(3,500,Phi=list(par=Phi,formula=~cov),p=p,design.data=ddi)
#    There are only 1000, because the last cohort is not included because it contains
#    no information
		xx$cov=cov$cov[1:1000]
		mod<-mark(xx,title="sim test",parameters=list(Phi=list(formula=~cov),p=list(formula=~1)),output=FALSE)
		results.real=rbind(results.real,mod$results$real$estimate)
		results.beta=rbind(results.beta,mod$results$beta$estimate)
	}
	return(list(real=results.real,beta=results.beta))
}

library(marked)
do.cjs.compare=function(Phi,p,num.cohorts,cohort.sizes,reps,...)
#
#  do.cjs.compare -  a simple example piece of code to show simulation/analysis loop
#            with a CJS model using crm and RMark
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
		xx=simcjs(num.cohorts,cohort.sizes,Phi=Phi,p=p)
		mod<-mark(xx,title="sim test",output=FALSE,...)
		results=rbind(results,mod$results$beta$estimate)
		results.se=rbind(results.se,mod$results$beta$se)
		mod<-crm(xx,hessian=TRUE)
		crmresults=rbind(crmresults,mod$beta)
		crmresults.se=rbind(crmresults.se, sqrt(diag(mod$vcv)))
		gc()
	}
	return(list(mark=list(beta=results,se=results.se),crm=list(beta=crmresults,se=crmresults.se)))
}
xx=do.cjs.compare(.9,.4,7,100,reps=25)
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}

