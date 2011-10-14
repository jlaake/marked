crm <- function(data,ddl=NULL,begin.time=1,model="cjs",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
 groups = NULL, time.intervals = NULL,debug=FALSE, method="nlminb", hessian=FALSE, accumulate=TRUE,chunk_size=1e7, 
 control=NULL,refit=1,itnmax=500,scale=NULL,autoscale=0,run=TRUE,...)
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
runmodel$model_data=NULL
if(runmodel$convergence!=0)cat("\n ******Model did not converge******\n")
if(runmodel$convergence==1)cat("\n Maximum number of iterations exceeded\n")
runmodel$model.parameters=model.parameters
cat(paste("\n Elapsed time in minutes: ",(proc.time()[3]-ptm[3])/60),"\n")
return(runmodel)
}
crm.wrapper <- function(model.list,data,ddl=NULL,...)
{
	for (i in 1:nrow(model.list))
	{
		model.parameters=list()
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
		model.name=paste(model.list[i,],collapse=".")
		cat("\n",model.name,"\n")
		mymodel=crm(data=data,ddl=ddl,model.parameters=model.parameters)
		assign(as.character(as.name(model.name)),mymodel)
		eval(parse(text=paste("save(",model.name,', file="',model.name,'.rda")',sep="")))
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
