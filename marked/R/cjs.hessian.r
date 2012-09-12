#' Compute variance-covariance matrix for fitted CJS model
#' 
#' A wrapper function that sets up call to hessian function to compute and then
#' invert the hessian.
#' 
#' @param model fitted CJS model from function crm
#' @param Phi.links vector of link function names for Phi parameters (currently unused)
#' @param p.links vector of link function names for p parameters (currently unused)
#' @export
#' @return variance-covariance matrix for specified model or the model
#' object with the stored vcv depending on whether the model has already been run
#' @author Jeff Laake <jeff.laake@@noaa.gov>
cjs.hessian=function(model,Phi.links=NULL, p.links=NULL)
{
	object=NULL
	if(!is.null(model$results))
	{
		object=model
		model=model$results
	}	
	scale=c(model$scale$phi,model$scale$p)
	assign(".markedfunc_eval", 0, envir = .GlobalEnv)
	vcv=hessian(cjs.lnl,model$beta*scale,model_data=model$model_data,Phi.links=NULL, p.links=NULL,all=FALSE)
	assign(".markedfunc_eval", 0, envir = .GlobalEnv)
	vcv=try(solvecov(vcv))
	if(class(vcv)[1]=="try-error")
	{
		warning("Unable to invert hessian")
		return(NULL)
	}
	vcv=vcv$inv/outer(scale,scale,"*")
	colnames(vcv)=names(model$beta)
	rownames(vcv)=names(model$beta)
	if(is.null(object))
	   return(vcv)
    else
	{
		model$beta.vcv=vcv
		object$results=model
		return(object)
	}
}
