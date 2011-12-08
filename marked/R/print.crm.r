#' Print model results or extract coefficients
#' 
#' Provides a printed simple summary of the model results or extracts the beta
#' coefficients from the model results.
#' 
#' @usage \method{coef}{crm}(object,...)
#'        \method{print}{crm}(x,...)
#' @S3method coef crm
#' @S3method print crm
#' @aliases print.crm coef.crm
#' @param x crm model result
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
   cat("\ncrm Model Summary\n")
   cat("\nNpar : ",length(x$beta))
   cat("\n-2lnL: ",x$neg2lnl)
   cat("\nAIC  : ",x$AIC)
   cat("\n\nBeta\n")
   print(coef(x))
   return(NULL)
}
coef.crm=function(object,...)
{
   beta=data.frame(Estimate=object$beta)
   if(!is.null(object$vcv))
   {
      beta$se=sqrt(diag(object$vcv))
      beta$lcl=beta$Estimate - 1.96*beta$se
      beta$ucl=beta$Estimate + 1.96*beta$se
   }
   return(beta)
}
