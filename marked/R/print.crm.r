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
