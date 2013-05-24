#' Print model results or extract coefficients
#' 
#' Provides a printed simple summary of the model results or extracts the beta
#' coefficients from the model results.
#' 
#' @usage \method{coef}{crm}(object,...)
#'        \method{print}{crm}(x,...)
#'        \method{print}{crmlist}(x,...)
#' @S3method coef crm
#' @S3method print crm
#' @S3method print crmlist
#' @aliases print.crm coef.crm print.crmlist
#' @param x crm model result or list of model results
#' @param object crm model result
#' @param ... generic arguments not used here
#' @return \code{print} prints a simple summary of the model to the screen and
#' returns NULL. \code{coef} returns a dataframe with estimates and standard
#' errors and confidence intervals if hessian=TRUE on model run.
#' @author Jeff Laake
#' @export 
#' @seealso \code{\link{crm}}
#' @keywords utility
print.crm=function(x,...)
{
   if(!is.null(x$results))x=x$results
   cat("\ncrm Model Summary\n")
   if(class(x)[2]=="mcmc")
	   cat("\nNpar : ",sum(sapply(x$beta,nrow)))   
   else
   {
	   cat("\nNpar : ",sum(sapply(x$beta,length)))
	   cat("\n-2lnL: ",x$neg2lnl)
	   cat("\nAIC  : ",x$AIC)
   }
   cat("\n\nBeta\n")
   print(coef(x))
   return(NULL)
}
coef.crm=function(object,...)
{
   if(class(object)[2]=="mcmc")
   {
	   beta=do.call(rbind,object$beta)
	   indices=grep("\\.",rownames(beta))
	   rownames(beta)[-indices]=paste(rownames(beta)[-indices],"(Intercept)",sep=".")
   }
   else
   {
	   beta=data.frame(Estimate=unlist(object$beta))
	   if(!is.null(object$beta.vcv))
	   {
		   beta$se=sqrt(diag(object$beta.vcv[1:length(beta$Estimate),1:length(beta$Estimate)]))
		   beta$lcl=beta$Estimate - 1.96*beta$se
		   beta$ucl=beta$Estimate + 1.96*beta$se
	   }
   }
   return(beta)
}
print.crmlist<-function(x,...)
{
	if(!is.null(x$model.table))
		print(x$model.table)
	else 
		cat("No model.table is available")
}
