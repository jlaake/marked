#' Automation of model runs
#' 
#' Some functions that help automate running a set of crm models based on parameter
#' specifications.
#' 
#' create.model.list creates all combinations of model specifications for the specified
#' set of parameters.  In the calling environment it looks for objects named parameter.xxxxxx where xxxxxx can
#' be anything. It creates a matrix with a column for each parameter and as many rows 
#' needed to create all combinations. This can be used as input to crm.wrapper.
#' 
#' crm.wrapper runs a sequence of crm models by constructing the call with the arguments
#' and the parameter specifications.  The parameter specifications can either be in the
#' local environment or in the environment of the named function models. The advantage of the
#' latter is that it is self-contained such that sets of parameter specifications can
#' be selected without possibility of being over-written or accidentally changed whereas 
#' with the former the set must be identified via a script and any in the environment will
#' be used which requires removing/recreating the set to be used.
#'  
#' @aliases  crm.wrapper create.model.list
#' @usage  crm.wrapper(model.list,data,ddl=NULL,models=NULL,base="",...)
#' 
#' create.model.list(parameters)
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param model.list matrix of model names contained in the environment of models function; each row is a model and each column is for a parameter and the value is the ime of first capture(release) occasion
#' @param models a function with a defined environment with model specifications as variables; values of model.list are some or all of those variables
#' @param base base value for model names
#' @param parameters character vector containing parameters in the model
#' @param ... aditional arguments passed to crm
#' @return create.model.list returns a matrix for crm.wrapper; crm.wrapper runs and stores models externally and has no return value
#' either "cjs" or "js" at present.
#' @author Jeff Laake
#' @export create.model.list
#' @export crm.wrapper
#' @seealso \code{\link{crm}}
#' @keywords models
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
		vec=ls(pattern=paste("^",n,"\\.",sep=""),envir=parent.frame())
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
