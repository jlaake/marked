#' Compute estimates of real parameters
#' 
#' Computes real estimates and their var-cov for a particular subset of 
#' parameters. The argument newdata may not work with all models. A better approach to 
#' compute real estimates for a subset of values or a new set of values is to specify a limited 
#' range of the values in ddl for each parameter. Make sure to include a complete set of values that spans
#' the factor levels and individual covariates used in the formulas for the model object or you will receive an
#' error that the number of columns in the design matrix does not match the number of beta parameters.  You cannot 
#' change the levels of any factor variable or modify the design data in anyway that changes the design matrix.
#' 
#' @usage \method{predict}{crm}(object,newdata=NULL,ddl=NULL,parameter=NULL,unique=TRUE,
#'                    vcv=FALSE,se=FALSE,chat=1,subset,select,...)
#' @param object model object;
#' @param newdata a dataframe for crm 
#' @param ddl list of dataframes for design data
#' @param parameter name of real parameter to be computed (eg "Phi")
#' @param unique TRUE if only unique values should be returned
#' @param vcv logical; if TRUE, computes and returns v-c matrix of real estimates
#' @param se logical; if TRUE, computes std errors and conf itervals of real estimates
#' @param chat over-dispersion value
#' @param subset logical expression using fields in real dataframe
#' @param select character vector of field names in real that you want to include 
#' @param merge default FALSE but if TRUE, the ddl for the parameter is merged (cbind) to the estimates
#' @param ... generic arguments not used here
#' @return A data frame (\code{real}) is returned if \code{vcv=FALSE};
#' otherwise, a list is returned also containing vcv.real: \item{real}{ data
#' frame containing estimates, and if vcv=TRUE it also contains
#' standard errors and confidence intervals} \item{vcv.real}{variance-covariance matrix of
#' real estimates}
#' @author Jeff Laake
#' @export
#' @examples
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,
#'    model.parameters=list(Phi=list(formula=~sex+time),p=list(formula=~1)),hessian=TRUE)
#' xx=predict(mod.Phisex.pdot,ddl=dipper.ddl)
#' xx
#' xx=predict(mod.Phisex.pdot,newdata=dipper[c(1,23),],vcv=TRUE)
#' xx
#' @keywords utility
predict.crm <-function(object,newdata=NULL,ddl=NULL,parameter=NULL,unique=TRUE,vcv=FALSE,se=FALSE,chat=1,subset=NULL,
                       select=NULL,real.ids=NULL,merge=FALSE,...)
{
  # if(object$model=="MVMSCJS")
  # {
  #   if(!is.null(newdata))message("\nargument newdata ignored for this model\n")
  #   if(!object$results$options$use.tmb)
  #     stop("Real predictions for MVMS model only works with TMB")
  #   if(is.null(ddl))
  #       stop("Must specify ddl argument")
  #   else  {
  #      emptyids=is.null(object$results$real.ids)&is.null(real.ids)
  #      if(!emptyids)
  #      {
  #        if(is.null(object$results$real.ids)&!is.null(real.ids))
  #          newids=TRUE
  #        else
  #          if(!is.null(object$results$real.ids)&is.null(real.ids))
  #            newids=TRUE
  #          else
  #            if(length(real.ids)!=length(object$results$real.ids))
  #              newids=TRUE
  #            else
  #              newids=!all(real.ids%in%object$results$real.ids)
  #      } else
  #        newids=FALSE
  #      if(is.null(object$results$reals) | newids)
  #         object=crm(object$data,ddl=ddl,model.parameters=object$model.parameters,optimize=FALSE,getreals=TRUE,
  #                          real.ids=real.ids,initial=object,use.tmb=TRUE,clean=FALSE,save.matrices=FALSE,check=FALSE)
  #      for(par in names(object$results$reals))
  #         object$results$reals[[par]]=cbind(ddl[[par]][ddl[[par]]$id %in% object$results$real.ids,],estimate=object$results$reals[[par]],se=object$results$reals.se[[par]])
  #      return(list(reals=object$results$reals))
  #   }
  # }
	if(!is.null(newdata))
	{
		if(is.data.frame(newdata))
		{
		  if(object$model=="CJS")
			   newdata$ch=paste(rep("1",object$data$nocc),collapse="")
		  else
		     newdata$ch=paste(rep(object$data$strata.labels[1],object$data$nocc),collapse="")
		  if(!is.null(object$data$group.covariates))
		  {
		    covs=apply(object$data$group.covariates,2,function(x) rep(as.character(x),each=nrow(newdata)))
		    if(nrow(object$data$group.covariates)>1)
		      for(i in 2:nrow(object$data$group.covariates))
		        newdata=rbind(newdata,newdata)
		    newdata=cbind(newdata,covs)
		  }
			newdata.proc=process.data(newdata,model=object$model,begin.time=object$data$begin.time,groups=names(object$data$group.covariates),strata.labels=object$data$strata.labels,accumulate=FALSE)
			ddl=make.design.data(newdata.proc,parameters=object$design.parameters)
			dml=create.dml(ddl,model.parameters=object$model.parameters,design.parameters=object$design.parameters)
		}else
			stop("Invalid newdata")
	} else
	{
		if(is.null(ddl))
		{
			if(!is.null(object$results$reals))
				return(object$results$reals)
			else{
				if(!is.null(object$results$model_data$ddl))
					ddl=object$results$model_data$ddl
				else
					stop("No ddl or real values available")
			}
		} else
		{ 
		  if(is.null(parameter))
			   dml=create.dml(ddl,model.parameters=object$model.parameters,design.parameters=ddl$design.parameters,chunk_size=1e7)  
		  else
		    dml=create.dml(ddl,model.parameters=object$model.parameters[parameter],design.parameters=ddl$design.parameters[parameter],chunk_size=1e7)   
		}
	}
	if(is.null(parameter))
	{
		results=NULL
		for (parameter in names(object$model.parameters))
			results[[parameter]]=compute_real(object,parameter,ddl,dml,unique,vcv,se,chat,subset=substitute(subset),select,include=object$model.parameters[[parameter]]$include,merge=merge)
		return(results)
	} else
		return(compute_real(object,parameter,ddl,dml,unique,vcv,se,chat,subset=substitute(subset),select,include=object$model.parameters[[parameter]]$include,merge=merge))	
}
