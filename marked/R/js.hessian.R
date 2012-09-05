#' Compute variance-covariance matrix for fitted JS model
#' 
#' A wrapper function that sets up call to hessian function to compute and then
#' invert the hessian.
#' 
#' @param model fitted JS model from function crm
#' @param nobstot number of unique caught at least once by group if applicable
#' @export
#' @return variance-covariance matrix for specified model
#' @author Jeff Laake <jeff.laake@@noaa.gov>
js.hessian=function(model,nobstot)
{
	scale=c(model$scale$phi,model$scale$p,model$scale$pent,model$scale$N)
	cat("Computing hessian\n")
	vcv=hessian(js.lnl,model$beta*scale,model_data=model$model_data,nobstot=nobstot)
	vcv=try(solvecov(vcv))
	if(class(vcv)[1]=="try-error")
	{
		warning("Unable to invert hessian")
		return(NULL)
	}
	vcv=vcv$inv/outer(scale,scale,"*")
	colnames(vcv)=names(model$beta)
	rownames(vcv)=names(model$beta)
	return(vcv)
}

