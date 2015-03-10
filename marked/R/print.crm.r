#' Print model results 
#' 
#' Provides a printed simple summary of the model results.
#' 
#' @usage \method{print}{crm}(x,...)
#' @param x crm model result or list of model results
#' @param ... generic arguments not used here
#' @return \code{print} prints a simple summary of the model to the screen and
#' returns NULL. 
#' @author Jeff Laake
#' @export 
#' @seealso \code{\link{crm}}
#' @keywords utility
print.crm=function(x,...)
{
   if(mode(x)=="character")x=load.model(x)
   if(!is.null(x$results))x=x$results
   if(class(x)[2]=="admb" & class(x)[3]=="cjs")
   {
	   class(x)[1]="admb"
       print(x)
   }
   else
   {
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
   }
   invisible(x)
}
#' Extract coefficients
#' 
#' Extracts the beta coefficients from the model results.
#' 
#' @usage \method{coef}{crm}(object,...)
#' @param object crm model result
#' @param ... generic arguments not used here
#' @return  returns a dataframe with estimates and standard
#' errors and confidence intervals if hessian=TRUE on model run.
#' @author Jeff Laake
#' @export 
#' @seealso \code{\link{crm}}
#' @keywords utility
coef.crm=function(object,...)
{
   if("results"%in%names(object)) object=object$results
   if(class(object)[2]=="mcmc")
   {
	   beta=do.call(rbind,object$beta)
	   indices=grep("\\.",rownames(beta))
	   rownames(beta)[-indices]=paste(rownames(beta)[-indices],"(Intercept)",sep=".")
   }
   else
   {
       if(class(object)[2]=="admb" & class(object)[3]=="cjs")
	   {
		   class(object)[1]="admb"
		   beta=coef(object)
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
   }
   return(beta)
}
#' Print model table from model list
#' 
#' @usage \method{print}{crmlist}(x,...)
#' @param x list of model results
#' @param ... generic arguments not used here
#' @return None#' @author Jeff Laake
#' @export coef.crm print.crm print.crmlist
#' @seealso \code{\link{crm}}
#' @keywords utility
print.crmlist<-function(x,...)
{
	if(!is.null(x$model.table))
		print(x$model.table)
	else 
		cat("No model.table is available")
}
