#' Simulates data from CJS model 
#'
#' Simulates capture-recapture data from CJS model of releases of multiple cohorts.
#' Allows general model for Phi and p or constant or time-specific model depending on input arguments Phi and p.
#' Optionally outputs file to outfile with .inp extension for input to program MARK.
#' 
#' The design.data dataframe contains a record for each animal for each time period it is contained in the study.  The default
#' fieldnames are id(unique # in cohort), cohort (cohort number from 1:num.cohorts), time (i:num.cohorts,
#' where i is cohort #) and age (assumed to increment by 1 each occasion and defaults with 0 initial age)
#' The function create.simdesign.data can be used to create an initial dataframe which can be supplemented 
#' with other covariates. If a formula is given but no design.data is not provided this function
#' calls create.simdesign.data to construct a default dataframe for the specified problem.
#'
#' For designation of models for Phi or p, there are 3 options here:
#'         constant model: Phi=list(par=value) 
#'         time model:     Phi=list(par=c(val1,val2,...valk)) where k=num.cohorts-1 is number of survival intervals
#'         general model:  Phi=list(par=c(val1,val2,...valk),formula=~yourformula)) k is number of cols in model matrix
#'
#' For first 2, identity link is assumed but logit link can be used. For formula, logit link is 
#' assumed and required.  The formula must match fields used in design.data.  Note that fields in 
#' design.data are numeric and the formula must specify as.factor to create a factor variable.  For
#' example a time model with 3 times and the same survivals in each can be specified as:
#'
#'         Phi=list(par=c(0.731,0.786)) or s=list(par=c(1,0.3),formula=~as.factor(time))
#' 
#' The above are not equivalent to Phi=list(par=c(1,0.3),formula=~time) which is 
#' where Phi=exp(1+0.3*time)/(1+exp(1+0.3*time)), Phi1=0.786, Phi2=0.831.  Note that survivals are
#' indexed by time at beginning of interval 1 for interal 1 to 2 and capture probabilities by
#' time of occasion.  Since these are truly recapture probabilities it starts with occasion 2, so
#' p=list(par=c(1,0.3),formula=~time) would give values: p2=0.832, p3=0.870
#
#' @param num.cohorts  number of cohorts; design is square with same number of c-r eventsas num.cohorts; number of recapture events is num.cohorts-1; it does not return ch values for last release cohort
#' @param cohort.sizes a scalar giving constant size of each cohort or a vector of sizes of length num.cohorts
#' @param Phi a list defining the survival model with the following elements (see details)
#'                    par      - a vector of parameter values
#'                    formula  - a formula to use with design.data to construct model
#'                    link     - link function used with model to create probabilites (not used at present)
#' @param p a list defining the capture probability model (same structure as Phi)
#' @param design.data a dataframe with design data that allows model construction for probabilities (see details).
#' @param outfile prefix name of the output file for the ch data. extension .inp is always added for MARK
#' @export simcjs
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @return dataframe with ch (capture history) as character
#' @examples 
#' library(RMark)
#' do.cjs <- function(Phi,p,reps,...)
#'#
#'#  do.cjs -  a simple example piece of code to show simulation/analysis loop
#'#            with a CJS model
#'#
#'#  Arguments:
#'#
#'#  Phi -  parameter list for Phi
#'#  p   -  parameter list for p
#'#  reps-  number of simulation reps
#'#  ... -  arguments passed to mark 
#'#
#'#
#'#  Value:
#'# 
#'#  results - for this simple example it will be a matrix of the real parameter estimates
#'#
#'#  Functions used: simcjs, mark
#'#
#'#
#'{
#'	results=NULL
#'	for(i in 1:reps)
#'	{
#'		cat("rep ",i,"\n")
#'		xx=simcjs(3,500,Phi=Phi,p=p)
#'		mod<-mark(xx,title="sim test",output=FALSE,...)
#'		results=rbind(results,mod$results$real$estimate)
#'	}
#'	return(results)
#'}
#' do.cjscov <-	function(Phi,p,reps)
#'#
#'#  do.cjscov -  a simple example piece of code to show simulation/analysis loop
#'#               with a CJS model including a single covariate
#'#
#'#  Arguments:
#'#
#'#  Phi -  parameter list for Phi
#'#  p   -  parameter list for p
#'#  reps-  number of simulation reps
#'#
#'#
#'#  Value:
#'# 
#'#  results - for this simple example it will be a list with 2 elements: a matrix of the real parameter estimates
#'#            and a matrix of beta estimates
#'#
#'#  Functions used: create.simdesign.data, simcjs, mark
#'#
#'#
#'{
#'	results.real=NULL
#'	results.beta=NULL
#'	for(i in 1:reps)
#'	{
#'		cat("rep ",i,"\n")
#'		cov=data.frame(id=1:1500,cov=rnorm(1500,0,1))
#'		ddi=create.simdesign.data(3,500,0)
#'		ddi=merge(cov,ddi,by="id",sort=FALSE)
#'		xx=simcjs(3,500,Phi=list(par=Phi,formula=~cov),p=p,design.data=ddi)
#' #    There are only 1000, because the last cohort is not included because it contains
#' #    no information
#'		xx$cov=cov$cov[1:1000]
#'		mod<-mark(xx,title="sim test",parameters=list(Phi=list(formula=~cov),p=list(formula=~1)),output=FALSE)
#'		results.real=rbind(results.real,mod$results$real$estimate)
#'		results.beta=rbind(results.beta,mod$results$beta$estimate)
#'	}
#'	return(list(real=results.real,beta=results.beta))
#'}
#'
#' library(marked)
#'do.cjs.compare=function(Phi,p,num.cohorts,cohort.sizes,reps,...)
#'#
#'#  do.cjs.compare -  a simple example piece of code to show simulation/analysis loop
#'#            with a CJS model using crm and RMark
#'#
#'#  Arguments:
#'#
#'#  Phi          parameter list for Phi
#'#  p            parameter list for p
#'#  num.cohorts  number of cohorts
#'#  cohort.sizes size or sizes of cohorts
#'#  reps         number of simulation reps
#'#  ...          additional optional arguments passed to mark
#'#
#'#
#'#  Value:
#'#
#'#  results - for this simple example it will be a matrix of the beta parameter estimates and
#'#	         a matrix of the standard errors
#'#
#'#  Functions used: simcjs, mark
#'#
#'#
#'{
#'	results=NULL
#'	results.se=NULL
#'	crmresults=NULL
#'	crmresults.se=NULL
#'	for(i in 1:reps)
#'	{
#'		cat("rep ",i,"\n")
#'		xx=simcjs(num.cohorts,cohort.sizes,Phi=Phi,p=p)
#'		mod<-mark(xx,title="sim test",output=FALSE,...)
#'		results=rbind(results,mod$results$beta$estimate)
#'		results.se=rbind(results.se,mod$results$beta$se)
#'		mod<-crm(xx,hessian=TRUE)
#'		crmresults=rbind(crmresults,mod$beta)
#'		crmresults.se=rbind(crmresults.se, sqrt(diag(mod$vcv)))
#'		gc()
#'	}
#'	return(list(mark=list(beta=results,se=results.se),crm=list(beta=crmresults,se=crmresults.se)))
#'}
#'xx=do.cjs.compare(.9,.4,7,100,reps=25)
#' 
simcjs <- function(num.cohorts=1,cohort.sizes,Phi=list(),p=list(),design.data=NULL,outfile=NULL)
{
#
# Setup default values for link and allow simple parameter vectors
#
	if(!is.list(Phi))
		Phi=list(par=Phi,link="identity")
	if(!is.list(p))
		p=list(par=p,link="identity")
	if(is.null(Phi$link))Phi$link="logit"
	if(is.null(p$link))p$link="logit"
#
# Check to make sure that Phi$par and p$par are defined
#
	if(is.null(Phi$par)) stop("Survival parameters are missing\n")
	if(is.null(p$par)) stop("Capture  parameters are missing\n")
#       
# Setup cohorts and cohort sizes
#
	ch=NULL
	if(num.cohorts>1 & length(cohort.sizes)==1)
		cohort.sizes=rep(cohort.sizes,num.cohorts)
	else
		if(missing(cohort.sizes))
			stop("\ncohort.sizes was not specified\n")
	    else
			if(num.cohorts !=length(cohort.sizes))
		        stop("\nNumber of cohorts doesn't match vector of cohort sizes\n")
#
# Check if p,s are using formula and then make sure design.data exists;
# if it doesn't exist create default design data
#
	if(!is.null(Phi$formula) | !is.null(p$formula))
		if(is.null(design.data))
			design.data=create.simdesign.data(num.cohorts,cohort.sizes,initial.age=0)
#
# Set up survival parameters
#
	if(is.null(Phi$formula))
	{
		sformula=FALSE
		if(length(Phi$par)==1)
			sconstant=TRUE
		else
		if(length(Phi$par)==(num.cohorts-1))
			sconstant=FALSE
		else
			stop("Incorrect number of survival parameters; doesn't match number of occasions-1\n")
		if(Phi$link=="identity")   
			Phi=Phi$par
		else
			Phi=exp(Phi$par)/(1+exp(Phi$par))
	}
	else
	{
		if(Phi$link!="logit")cat("Note: Logit link will be used with formula\n")
		sformula=TRUE
		smat=model.matrix(Phi$formula,design.data[design.data$cohort<num.cohorts&design.data$time<num.cohorts,])
		if(dim(smat)[2]!=length(Phi$par))
			stop(paste("Dimension of model matrix",dim(smat)[2],"is not consistent with length of survival parameter vector",length(Phi$par),"\n"))
	}
#
# Set up capture parameters
#
	if(is.null(p$formula))
	{
		pformula=FALSE
		if(length(p$par)==1)
			pconstant=TRUE
		else
		if(length(p$par)==(num.cohorts-1))
			pconstant=FALSE
		else
			stop("Incorrect number of capture parameters; doesn't match number of occasions-1\n")
		if(p$link=="identity")   
			p=p$par
		else
			p=exp(p$par)/(1+exp(p$par))
	}
	else
	{
		if(p$link!="logit")cat("Note: Logit link will be used with formula\n")
		pformula=TRUE
		pmat=model.matrix(p$formula,design.data[design.data$cohort<num.cohorts&design.data$time>design.data$cohort,])
		if(dim(pmat)[2]!=length(p$par))
			stop(paste("Dimension of model matrix",dim(pmat)[2],"is not consistent with length of capture parameter vector",length(p$par),"\n"))
	}
#
# Loop over cohorts and generate simulated data for each cohort
#
	for (i in 1:(num.cohorts-1))
	{
#
#  Setup cohort specific survival parameters
#
		if(!sformula)
		{
			if(!sconstant)  
				svalues=Phi[i:(num.cohorts-1)]
			else
				svalues=Phi
		} 
		else
		{
			svalues=smat[design.data$cohort[design.data$time<num.cohorts]==i,] %*% Phi$par
			svalues=exp(svalues)/(1+exp(svalues))
			select=design.data$cohort==i&design.data$time>=i&design.data$time<num.cohorts
			svalues=tapply(as.vector(svalues),list(as.factor(design.data$id[select]),as.factor(design.data$time[select])),sum)
			if(dim(svalues)[1]!=cohort.sizes[i]) 
				stop("Number of rows in Phi model matrix not consistent with number in cohort\n")
		}
#
#  Setup cohort specific capture parameters
#
		if(!pformula)
		{
			if(!pconstant)
				pvalues=p[i:(num.cohorts-1)]       
			else
				pvalues=p
		}
		else
		{
			pvalues=pmat[design.data$cohort[design.data$time>design.data$cohort&design.data$cohort<num.cohorts]==i,] %*% p$par
			pvalues=exp(pvalues)/(1+exp(pvalues))
			select=design.data$cohort==i&design.data$time>i
			pvalues=tapply(as.vector(pvalues),list(as.factor(design.data$id[select]),as.factor(design.data$time[select])),sum)
			if(dim(pvalues)[1]!=cohort.sizes[i]) 
				stop("Number of rows in p model matrix not consistent with number in cohort\n")
		}
#
#  Create simulated capture histories for the cohort and append to ch; if this is 
#  not the first cohort, the appropriate number of 0's are used as prefix to the
#  cohort capture history.
#
		ch=c(ch,paste(rep(paste(rep("0",i-1),collapse=""),cohort.sizes[i]),simcjs.cohort(cohort.sizes[i],svalues,pvalues,num.cohorts-i+1),sep=""))   
	}
#
#  If a filename specified output ch to file for input to MARK
#
	if(!is.null(outfile))
		write.table(paste(ch,rep(" 1 ;",length(ch))),paste(outfile,".inp",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
#
#  Return ch vector as value 
#
	return(data=data.frame(ch=I(ch)))
}
simcjs.cohort <-
		function(N,Phi,p,nocc)
#
# simcjs.cohort - simulates capture histories for a single cohort of N release animals 
#         for nocc occasions using the CJS model 
#
# Arguments: 
#
#  N        - size of cohort
#  Phi      - survival probability (can be either scalar(constant model), 
#             vector (time model - length = nocc-1 or heterogeneity model - length=N)
#             or matrix (N by nocc-1)
#  p        - capture probability (possibilities same as s)
#  nocc     - number of occasions (ch length is nocc)
#
# Value:
#  ch       - character vector of capture histories 
#
# Functions called:  create.parmat
#
{
#
# Setup parameter matrix depending on what was input; it can be a scalar (constant model), a vector with
# nocc-1 values (time model) or N values (heterogeneity) or a matrix (N rows by nocc-1 cols) for a completely 
# general specification.  In each case it is transformed to a matrix.
#
	Phi=create.parmat(Phi,nocc,N)
	p=create.parmat(p,nocc,N)
#
# alive is matrix of who is alive (row=id, col=occasion) 
#
   alive=matrix(rbinom(N*(nocc-1),1,Phi),nrow=N,ncol=nocc-1)
   if(nocc>2) alive=t(apply(alive,1,cumprod))
#
#  Generate Bernoulli rv for seen
#
   seen=matrix(rbinom(N*(nocc-1),1,p),nrow=N,ncol=nocc-1)
#
#  Create capture history 
#
   ch=apply(cbind(rep("1",N),alive*seen),1,paste,collapse="")
#
#    Return a vector of capture histories with a 1 ; appended for input to MARK;
#    Vector can be written to file with command
#    write.table(ch,filename,row.names=FALSE,col.names=FALSE,quote=FALSE)
#
	return(ch)
}


