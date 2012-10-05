#' Compute estimates of real parameters
#' 
#' Computes real estimates and their var-cov for a particular subset of 
#' parameters.
#' 
#' @param object model object
#' @param newdata a dataframe for crm 
#' @param parameter name of real parameter to be computed (eg "Phi" or "p")
#' @param unique TRUE if only unique values should be returned
#' @param vcv logical; if TRUE, computes and returns v-c matrix of real estimates
#' @param se logical; if TRUE, computes std errors and conf itervals of real estimates
#' @param chat over-dispersion value
#' @param subset logical expression using fields in real dataframe
#' @param select character vector of field names in real that you want to include 
#' @return A data frame (\code{real}) is returned if \code{vcv=FALSE};
#' otherwise, a list is returned also containing vcv.real: \item{real}{ data
#' frame containing estimates, and if vcv=TRUE it also contains
#' standard errors and confidence intervals} \item{vcv.real}{variance-covariance matrix of
#' real estimates}
#' @author Jeff Laake
#' @S3method predict crm
#' @export
#' @examples
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~sex+time),p=list(formula=~1)),hessian=T)
#' xx=predict(mod.Phisex.pdot)
#' xx
#' xx=predict(mod.Phisex.pdot,newdata=dipper[c(1,23),],vcv=TRUE)
#' xx
#' @keywords utility
predict.crm <-function(object,newdata=NULL,ddl=NULL,parameter=NULL,unique=TRUE,vcv=FALSE,se=FALSE,chat=1,subset,select)
{
	if(!is.null(newdata))
	{
		if(is.data.frame(newdata))
		{
			newdata$ch=paste(rep("1",object$data$nocc),collapse="")
			newdata.proc=process.data(newdata,model=object$model,begin.time=object$data$begin.time,groups=names(object$data$group.covariates),accumulate=FALSE)
			ddl=make.design.data(newdata.proc,parameters=object$design.parameters)
			dml=create.dml(ddl,model.parameters=object$model.parameters,design.parameters=object$design.parameters)
		    object$results$model_data$Phi.dm=dml$Phi
			object$results$model_data$p.dm=dml$p		
		}else
			stop("Invalid newdata")
	}
	if(is.null(ddl))
		return(object$results$reals)
	else
	{
		if(is.null(parameter) )
		{
			results=NULL
			for (parameter in names(object$model.parameters))
				results[[parameter]]=compute.real(object,parameter,ddl,unique,vcv,se,chat,subset=substitute(subset),select)
			return(results)
		} else
			return(compute.real(object,parameter,ddl,unique,vcv,se,chat,subset=substitute(subset),select))		
	}
}