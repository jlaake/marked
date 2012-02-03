#' Capture-recapture model fitting function
#' 
#' Fits user specified models to some types of capture-recapture wholly in R
#' and not with MARK.  A single function that processes data, creates the
#' design data, makes the crm model and runs it
#' 
#' This function is operationally similar to the function \code{mark in RMark}
#' in that is is a shell that calls several other functions to perform the
#' following steps: 1) \code{\link{process.data}} to setup data and parameters
#' and package them into a list (processed data),2)
#' \code{\link{make.design.data}} to create the design data for each parameter
#' in the specified model, 3) \code{\link{create.dm}} to create the design
#' matrices for each parameter based on the formula provided for each
#' parameter, 4) call to the specific function for model fitting (now either
#' \code{\link{cjs}} or \code{\link{js}}). As with \code{mark} the calling
#' arguments for \code{crm} are a compilation of the calling arguments for each
#' of the functions it calls (with some arguments renamed to avoid conflicts).
#' If data is a processed dataframe (see \code{\link{process.data}}) then it
#' expects to find a value for \code{ddl}.  Likewise, if the data have not been
#' processed, then \code{ddl} should be NULL.  This dual calling structure
#' allows either a single call approach for each model or alternatively for the
#' data to be processed and the design data (\code{ddl}) to be created once and
#' then a whole series of models can be analyzed without repeating those steps.
#' 
#' There are some optional arguments that can be used to set initial values and
#' control other aspects of the optimization.  The optimization is done with
#' the R function \code{optim} and the arguments \code{method} and
#' \code{hessian} are described with the help for that function.  In addition,
#' any arguments not matching those for \code{cjs} (the ...) are passed to
#' \code{optim} allowing any of the other parameters to be set.  If you set
#' \code{debug=TRUE}, then at each function evaluation (\code{\link{cjs.lnl}}
#' the current values of the parameters and -2*log-likelihood value are output.
#' 
#' In the current implementation, a logit link is used to constrain the
#' parameters in the unit interval (0,1) except for probability of entry which
#' uses an mlogit and N which uses a log link.  These could be generalized to
#' use other link functions. Following the notation of MARK, the parameters in
#' the link space are referred to as \code{beta} and those in the actual
#' parameter space of \code{Phi} and \code{p} as reals.
#' 
#' Initial values can be set in 2 ways.  To set a baseline intial value for the
#' intercept of \code{Phi} \code{p} set those arguments to some real value in
#' the open interval (0,1). All non-intercept beta parameters are set to zero.
#' Alternatively, you can specify in \code{initial}, a vector of initial values
#' for the beta parameters (on the logit scale).  This is most easily done by
#' passing the results from a previous model run using the result list element
#' \code{beta} as described below.  The code will match the names of the
#' current design matrix to the names in \code{beta} and use the appropriate
#' initial values. Any non-specified values are set to 0.  If there are no
#' names associated with the \code{initial} vector then they are simply used in
#' the specified order. If you do not specify initial values it is equivalent
#' to setting \code{Phi} and \code{p} to 0.5.
#' 
#' If you have a study with unequal time intervals between capture occasions,
#' then these can be specified with the argument \code{time.intervals}.
#' 
#' The argument \code{accumulate} defaults to \code{TRUE}.  When it is
#' \code{TRUE} it will accumulate common capture histories that also have
#' common design and common fixed values (see below) for the parameters.  This
#' will speed up the analysis because in the calculation of the likelihood
#' (\code{\link{cjs.lnl}} it loops over the unique values. In general the
#' default will be the best unless you have many capture histories and are
#' using many individual covariate(s) in the formula that would make each entry
#' unique. In that case there will be no effect of accumulation but the code
#' will still try to accumulate. In that particular case by setting
#' \code{accumulate=FALSE} you can skip the code run for accumulation.
#' 
#' Most of the arguments controlling the fitted model are contained in lists in
#' the arguments \code{model.parameters} and \code{design.parameters} which are
#' similar to their counterparts in \code{mark inb RMark}. Each is a named list
#' with the names being the parameters in the model (e.g., Phi and p in "cjs"
#' and "Phi","p","pent","N" in "js"). Each named element is also a list
#' containing various values defining the design data and model for the
#' parameter. The elements of \code{model.parameters} can include
#' \code{formula} which is an R formula to create the design matrix for the
#' parameter and \code{fixed} is a matrix of fixed values as described below.
#' The elements of \code{design.parameters} can include \code{time.varying},
#' \code{fields}, \code{time.bins},\code{age.bins}, and \code{cohort.bins}. See
#' \code{\link{create.dmdf}} for a description of the first 2 and
#' \code{\link{create.dm}} for a description of the last 3.
#' 
#' Real parameters can be set to fixed values using \code{fixed=x} where x is a
#' matrix with 3 columns and any number of rows.  The first column specifies
#' the particular animal (capture history) as the row number in the dataframe
#' \code{x}.  The second specifies the capture occasion number for the real
#' parameter to be fixed.  For \code{Phi} and \code{pent} these are 1 to
#' \code{nocc}-1 and for \code{p} they are 2 to \code{nocc} for "cjs" and 1 to
#' \code{nocc} for "js". This difference is due to the parameter labeling by
#' the beginning of the interval for Phi (e.g., survival from occasion 1 to 2)
#' and by the occasion for \code{p}.  For "cjs" \code{p} is not estimated for
#' occasion 1. The third element in the row is the real value in the closed
#' unit interval [0,1] for the fixed parameter.  This approach is completely
#' general allowing you to fix a particular real parameter for a specific
#' animal and occasion but it is a bit kludgy and later some other
#' functionality will be included.  To set all of the real values for a
#' particular occasion you can use the following example with the dipper data
#' as a template:
#' 
#' \code{model.parameters=list(Phi=list(,fixed=cbind(1:nrow(dipper),rep(2,nrow(dipper)),rep(1,nrow(dipper)))))}
#' 
#' The above sets \code{Phi} to 1 for the interval between occasions 2 and 3
#' for all animals. At present there is no modification of the parameter count
#' to address fixing of real parameters.
#' 
#' Be cautious with this function at present.  It does not include many checks
#' to make sure values are in the specified range of the data.  Normally this
#' would not be a big problem but because \code{\link{cjs.lnl}} calls an
#' external FORTRAN subroutine, if it gets a subscript out of bounds, it will
#' cause R to terminate.  So make sure to save your workspace frequently if you
#' use this function in its current implementation.
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model ("cjs" or "js" at present)
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional vector of initial values for beta parameters; if
#' named from previous analysis only relevant values are used
#' @param groups Vector of names factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param method optimization method for function \code{optim}
#' @param debug if TRUE, shows optimization output
#' @param hessian if TRUE, computes v-c matrix using hessian
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce
#' computation
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories and design matrices; amount used is 8*chunk_size/1e6 MB (default
#' 80MB)
#' @param control control string for optimization functions
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param scale vector of scale values for parameters
#' @param ... optional arguments passed to js or cjs and optimx
#' @return crm model object with class=("crm",submodel) where submodel is
#' either "cjs" or "js" at present.
#' @author Jeff Laake
#' @export crm
#' @export create.model.list
#' @export crm.wrapper
#' @useDynLib marked
#' @seealso \code{\link{cjs}}, \code{\link{js}},
#' \code{\link{make.design.data}},\code{\link{process.data}}
#' @keywords models
#' @examples
#' 
#' # cormack-jolly-seber model
#' # fit 3 cjs models with crm
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phit.pt=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)))
#' mod.Phit.pt
#' mod.Phidot.pdot=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~1),p=list(formula=~1)))
#' mod.Phidot.pdot
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)))
#' mod.Phisex.pdot
#' # fit same 3 models with calls to mark
#' #require(RMark)
#' #mod0=mark(dipper,model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)),output=FALSE)
#' #summary(mod0,brief=TRUE)
#' #mod1=mark(dipper,model.parameters=list(Phi=list(formula=~1),p=list(formula=~1)),output=FALSE)
#' #summary(mod1,brief=TRUE)
#' #mod2<-mark(dipper,groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)),output=FALSE)
#' #summary(mod2,brief=TRUE)
#' # jolly seber model
#' crm(dipper,model="js",groups="sex",accumulate=FALSE)
#' #require(RMark)
#' #mark(dipper,model="POPAN",groups="sex")
#' 
crm <- function(data,ddl=NULL,begin.time=1,model="cjs",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
 groups = NULL, time.intervals = NULL,debug=FALSE, method="nlminb", hessian=FALSE, accumulate=TRUE,chunk_size=1e7, 
 control=NULL,refit=1,itnmax=5000,scale=NULL,autoscale=0,run=TRUE,...)
{
# -----------------------------------------------------------------------------------------------------------------------
# crm -  a single function that processes data, creates the design data, makes the crm model and runs it.
#
# Arguments:
#
#  data                 - either the raw data which is a dataframe with at least one column named 
#                         ch which is a character field containing the capture history or a processed dataframe
#  ddl                  - design data list which contains an element for each parameter type; if NULL it is created
#  begin.time           - time of first capture(release) occasion
#  model                - type of c-r model (eg cjs, js at present) 
#  title                - a title for the analysis 
#  model.parameters     - list of parameter model specifications
#  design.parameters    - specification of any grouping variables for design data for each parameter
#  initial              - vector of initial values for beta parameters
#  groups               - list of factors for creating groups
#  time.intervals       - intervals of time between the capture occasions
#  debug                - if TRUE, shows optimization output
#  chunk_size           - specifies amount of memory to use in accumulating capture histories
#                             use is 8*chunk_size/1e6 MB (default 80MB)
#  refit                - number of times to refit the model if it doesn't converge
#  ...                  - additional arguments passed to cjs and js
#
#  Value: 
#
#  model - a crm object model containing output and extracted results
#
#  Functions used: process.data, create.dmdf, create.dm, cjs, js 
# 
# -----------------------------------------------------------------------------------------------------------------------
ptm=proc.time()	
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
# 
if(is.null(data$data))
{
   if(!is.null(ddl))
   {
      cat("Warning: specification of ddl ignored, as data have not been processed\n")
      ddl=NULL
   }
   cat("\n Processing data\n")
   flush.console()
   data.proc=process.data(data,begin.time=begin.time, model=model,mixtures=1, 
                          groups = groups, age.var = NULL, initial.ages = NULL, 
                          age.unit = NULL, time.intervals = time.intervals,nocc=NULL)
}   
else
   data.proc=data
#
# If the design data have not been constructed, do so now
#
if(is.null(ddl)) 
{
	cat("\n Creating design data. This can take awhile.\n")
	flush.console()
	ddl=make.design.data(data.proc,design.parameters)
}
#
# Setup parameter list
#
number.of.groups=1
if(!is.null(data.proc$group.covariates))number.of.groups=nrow(data.proc$group.covariates)
par.list=setup.parameters(data.proc$model,check=TRUE)
parameters=setup.parameters(data.proc$model,model.parameters,data$nocc,number.of.groups=number.of.groups)
parameters=parameters[par.list]
#
# Create design matrices for each parameter
#
dml=vector("list",length=length(parameters))
names(dml)=names(parameters)
for (i in 1:length(parameters))
{
   pn=names(parameters)[i]
   if(is.null(parameters[[i]]$formula)) parameters[[i]]$formula=~1
   cat("\n Creating design matrix for parameter ",pn," with formula ",paste(as.character(parameters[[i]]$formula),sep=""),"\n")
   flush.console()
   dml[[i]]=create.dm(ddl[[pn]],parameters[[i]]$formula,design.parameters[[pn]]$time.bins,
                                design.parameters[[pn]]$cohort.bins,design.parameters[[pn]]$age.bins,chunk_size=chunk_size,parameters[[i]]$remove.intercept)
}
if(!run) return(dml)
#
# Call estimation function
#
if(autoscale==0)
{
	control$eval.max=itnmax
    if(model=="cjs")
       runmodel=cjs(data.proc,ddl,dml,model_data=NULL,parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
		          refit=refit,control=control,itnmax=itnmax,scale=scale,...)
    else
       runmodel=js(data.proc,ddl,dml,parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
		          refit=refit,control=control,itnmax=itnmax,scale=scale,...)
}else
{
	cat("\n Run to compute scale:")
	scale=1
	control$eval.max=autoscale
	if(model=="cjs")
		runmodel=cjs(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=FALSE,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=0,control=control,itnmax=autoscale,scale=scale,...)
	else
		runmodel=js(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=FALSE,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=0,control=control,itnmax=autoscale,scale=scale,...)
	scale=abs(1/runmodel$beta)
	initial=runmodel$beta
	control$eval.max=itnmax
	cat("\n\n Fitting model:")
	if(model=="cjs")
		runmodel=cjs(data.proc,ddl,dml,model_data=runmodel$model_data,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=refit,control=control,itnmax=itnmax,scale=scale,...)
	else
		runmodel=js(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=refit,control=control,itnmax=itnmax,scale=scale,...)
}
#
# Return fitted MARK model object or if external, return character string with same class and save file
#
if(runmodel$convergence!=0)
{
	cat("\n ******Model did not converge******\n")
    msg=attr(runmodel$mod,"details")[[1]]$message
	if(is.null(msg)) msg="Exceeded maximum number of iterations"
    cat(paste(msg,"\n"))
}
runmodel$model.parameters=model.parameters
cat(paste("\n Elapsed time in minutes: ",(proc.time()[3]-ptm[3])/60),"\n")
return(runmodel)
}
crm.wrapper <- function(model.list,data,ddl=NULL,models=NULL,base="",...)
{
	for (i in 1:nrow(model.list))
	{
		model.parameters=list()		
		if(is.null(models))
		{
			for(j in 1:ncol(model.list))
			{
				if(!is.list(eval(parse(text=model.list[i,j]),envir=parent.frame())[[1]]))
					model.parameters[[names(model.list)[j]]]=eval(parse(text=(as.character(model.list[i,j]))),envir=parent.frame())
			}
			for(j in 1:ncol(model.list))
			{
				if(is.list(eval(parse(text=model.list[i,j]),envir=parent.frame())[[1]]))
					model.parameters=c(model.parameters,eval(parse(text=(as.character(model.list[i,j]))),envir=parent.frame()))
			}	
		} else
		{
			model.parameters=models(model.list[i,])$model.parameters
		}
		model.name=paste(model.list[i,],collapse=".")
		cat("\n",model.name,"\n")
		mymodel=crm(data=data,ddl=ddl,model.parameters=model.parameters,...)
		assign(as.character(as.name(model.name)),mymodel)
		eval(parse(text=paste("save(",model.name,', file="',base,model.name,'.rda")',sep="")))
	}	
}
create.model.list<-function(parameters)
{
	model.list=list()
	for(n in parameters)
	{
		vec=ls(pat=paste("^",n,"\\.",sep=""),envir=parent.frame())
		if(length(vec)>0)
			for (i in 1:length(vec))
			{
				if(eval(parse(text=paste("is.list(",vec[i],")",sep="")),envir=parent.frame()))
				{
					if(eval(parse(text=paste("!is.null(",vec[i],"$formula)",sep="")),envir=parent.frame()) |
							eval(parse(text=paste("!is.null(",vec[i],"[[1]]$formula)",sep="")),envir=parent.frame()))
						model.list[[n]]=c(model.list[[n]],vec[i])
				}
			}
		
	}
	if(length(model.list)==0)
		stop("\nNo model specifications found. Use parameter.description notation (e.g., Phi.time)\n")
	if(length(model.list)>1)
	{
		model.list=expand.grid(model.list)
		for (j in 1:dim(model.list)[2])
			model.list[,j]=as.character(model.list[,j])
	}
	else
		model.list=as.data.frame(model.list)
	return(model.list)
}
#
#
# solvecov code was taken from package fpc: Christian
# Hennig chrish@@stats.ucl.ac.uk http://www.homepages.ucl.ac.uk/~ucakche/
solvecov=function (m, cmax = 1e+10)
# from package fpc
{
	options(show.error.messages = FALSE)
	covinv <- try(solve(m))
	if (class(covinv) != "try-error")
		coll = FALSE
	else {
		p <- nrow(m)
		cove <- eigen(m, symmetric = TRUE)
		coll <- TRUE
		if (min(cove$values) < 1/cmax) {
			covewi <- diag(p)
			for (i in 1:p) if (cove$values[i] < 1/cmax)
					covewi[i, i] <- cmax
				else covewi[i, i] <- 1/cove$values[i]
		}
		else covewi <- diag(1/cove$values, nrow = length(cove$values))
		covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
	}
	options(show.error.messages = TRUE)
	out <- list(inv = covinv, coll = coll)
}

