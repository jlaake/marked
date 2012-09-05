#' Compute estimates of real parameters
#' 
#' Computes real estimates and var-cov from design matrix (design) and
#' coefficients (beta) using specified link functions
#' 
#' The estimated real parameters can be derived from the estimated beta
#' parameters, a completed design matrix, and the link function specifications.
#' MARK produces estimates of the real parameters, se and confidence intervals
#' but there are at least 2 situations in which it is useful to be able to
#' compute them after running the analysis in MARK: 1) adjusting confidence
#' intervals for estimated over-dispersion, and 2) making estimates for
#' specific values of covariates.  The first case is done in
#' \code{\link{get.real}} with a call to this function.  It is done by
#' adjusting the estimated standard error of the beta parameters by multiplying
#' it by the square root of \code{chat} to adjust for over-dispersion.  A
#' normal 95% confidence interval is computed for the link estimate (estimate
#' +/- 1.96*se) and this is then back-transformed to the real parameters using
#' \code{\link{inverse.link}} with the appropriate inverse link function for
#' the parameter to construct a 95% confidence interval for the real parameter.
#' There is one exception. For parameters using the \code{mlogit}
#' transformation, a \code{logit} transformation of each individual real Psi
#' and its se are used to derive the confidence interval. The estimated
#' standard error for the real parameter is also scaled by the square root of
#' the over-dispersion constant \code{chat} stored in \code{model$chat}. But,
#' the code actually computes the variance-covariance matrix rather than
#' relying on the values from the MARK output because real estimates will
#' depend on any individual covariate values used in the model which is the
#' second reason for this function.
#' 
#' New values of the real parameter estimates can easily be computed by simply
#' changing the values of the covariate values in the design matrix and
#' computing the inverse-link function using the beta parameter estimates.  The
#' covariate values to be used can be specified in one of 2 ways. 1) Prior to
#' making a call to this function, use the functions
#' \code{\link{find.covariates}} to extract the rows of the design matrix with
#' covariate values and either fill in those values aautomatically with the
#' options provided by \code{\link{find.covariates}} or edit those values to be
#' the ones you want and then use \code{\link{fill.covariates}} to replace the
#' values into the design matrix and use it as the value for the argument
#' \code{design}, or 2) automate this step by specifying a value for the
#' argument \code{data} which is used to take averages of the covariate values
#' to fill in the covariate entries of the design matrix.  In computing real
#' parameter estimates from individual covariate values it is important to
#' consider the scale of the individual covariates. By default, an analysis
#' with MARK will standardize covariates by subtracting the mean and dividing
#' by the standard deviation of the covariate value. However, in the
#' \code{RMark} library all calls to MARK.EXE do not standardize the covariates
#' and request real parameter estimates based on the mean covariate values.
#' This was done because there are many instances in which it is not wise to
#' use the standardization implemented in MARK and it is easy to perform any
#' standardization of the covariates with R commands prior to fitting the
#' models.  Also, with pre-standardized covariates there is no confusion in
#' specifying covariate values for computation of real estimates.  If the model
#' contains covariates and the argument \code{design} is not specified, the
#' design matrix is extracted from \code{model} and all individual covariate
#' values are assigned their mean value to be consistent with the default in
#' the MARK analysis.
#' 
#' If a value for \code{beta} is given, those values are used in place of the
#' estimates \code{model$results$beta$estimate}.
#' 
#' @param model MARK model object
#' @param parameter name of real parameter to be computed (eg "Phi" or "p")
#' @param unique TRUE if only unique values should be returned
#' @param vcv logical; if TRUE, computes and returns v-c matrix of real estimates
#' @param chat over-dispersion value
#' @param subset logical expression using fields in real dataframe
#' @param select character vector of field names in real that you want to include 
#' @return A data frame (\code{real}) is returned if \code{vcv=FALSE};
#' otherwise, a list is returned also containing vcv.real: \item{real}{ data
#' frame containing estimates, and if vcv=TRUE it also contains
#' standard errors and confidence intervals} \item{vcv.real}{variance-covariance matrix of
#' real estimates}
#' @author Jeff Laake
#' @export
#' @keywords utility
compute.real <-function(model,parameter,unique=FALSE,vcv=FALSE,chat=1,subset,select)
{
  links=c(Phi="logit",p="logit",pent="mlogit",N="log")
  link=links[parameter]
  if(!missing(subset))
  {
	  sub = substitute(subset)
	  indices = eval(sub, model$real)
  }
  else
	  indices = TRUE
  if(!missing(subset) | !missing(select))
  {
	  if(missing(select))
	     df=model$real[indices,,drop=FALSE]
      else
         df=model$real[indices,select,drop=FALSE]
   }else
	  df=model$real
# get dm
	design=model$model_data[[paste(parameter,"dm",sep=".")]][model$model_data$dmrec[model$real$rec],,drop=FALSE]
	design=design[indices,,drop=FALSE]
# if unique, find unique dm values
  if(unique)
  {
	  indices=which(!duplicated(design))
	  df=df[indices,,drop=FALSE]
	  design=design[indices,,drop=FALSE]
  }
#
# Check to make sure beta and DM match
#
beta=model$beta[grep(paste(parameter,":",sep=""),names(model$beta))]
if(dim(design)[2]!=length(beta))
   stop("Mismatch between number of design columns and length of beta")
#
# Set indices for real parameters that have been fixed and at any boundaries
#
#if(!is.null(model$fixed))
#{
#   if(is.null(model$simplify))
#      fixedparms=(1:dim(design)[1])%in%model$fixed$index
#   else
#      fixedparms=(1:dim(design)[1])%in%model$simplify$pim.translation[model$fixed$index]
#}
#else
   fixedparms=rep(FALSE,dim(design)[1])
#boundaryparms=model$results$real$se==0 & !fixedparms
fixedvalues=rep(NA,nrow(design))
#fixedvalues[model$simplify$pim.translation[model$fixed$index]]=model$fixed$value
#
#  Compute real parameters; if neither se or vcv then return vector of real parameters
#
   real=convert.link.to.real(design%*%beta,links=link,fixed=fixedvalues)
#
#  Set fixed real parameters to their fixed values
#
   real[fixedparms]=fixedvalues[fixedparms]
#  If no vcv requested, return result
if(!vcv)return(cbind(df,data.frame(estimate=real)))
if(is.null(model$vcv))stop("vcv matrix not available in model")
#
# Compute vc matrix for real parameters which needs to be modified for
# any mlogit parameters
#
indices=grep(paste(parameter,":",sep=""),colnames(model$vcv))
if(link!="mlogit")
{
	deriv.real=deriv_inverse.link(real,design,link)
	vcv.real=deriv.real%*%model$vcv[indices,indices]%*%t(deriv.real)
} else
#
# If vcv=TRUE, compute v-c matrix and std errors of real estimates
# To handle any mlogit parameters compute pseudo-real estimates using log in place of mlogit
#
{
	pseudo.real=as.vector(convert.link.to.real(design%*%beta,links="log"))
    pseudo.real[fixedparms]=fixedvalues[fixedparms]
    pseudo.real[fixedparms]=exp(pseudo.real[fixedparms])
#
#   Compute first derivatives of pseudo-real (for any mlogit parameter)
#   estimates with respect to beta parameters
#
    deriv.pseudo=deriv_inverse.link(pseudo.real,design,"log")
    deriv.pseudo[fixedparms,]=0
    vcv.pseudo=deriv.pseudo%*%model$vcv[indices,indices]%*%t(deriv.pseudo)
#
#    Apply chain rule to get variance of real parameters which has mlogits
#    expressed as zi/(1+z1+...zk) where k is number of mlogit components-1 and
#    non-mlogits are expressed as zi.
#    bottom is either 1 for non-mlogits and the sum for mlogits
#    pbottom is partial with respect to zi
#
    links=rep("mlogit(1)",length(pseudo.real))
    mlogits=outer(links,links,function(x,y)as.numeric(x==y))*as.numeric(substr(links,1,6)=="mlogit")
    pbottom=matrix(0,nrow=dim(vcv.pseudo)[1],ncol=dim(vcv.pseudo)[1]) + mlogits
    bottom=diag(nrow=nrow(vcv.pseudo))*(1-as.numeric(substr(links,1,6)=="mlogit"))+
       mlogits + pbottom*apply(pbottom*pseudo.real,2,sum)
    deriv.pseudo=(diag(nrow=dim(vcv.pseudo)[1])*bottom-pseudo.real*pbottom)/bottom^2
    deriv.pseudo[is.nan(deriv.pseudo)]=0
    vcv.real=deriv.pseudo%*%vcv.pseudo%*%t(deriv.pseudo)
}
vcv.real=chat*vcv.real
#
# Compute conf interval taking into account use of logit transform for mlogit
# and any 0-1 link (loglog,cloglog,sin,logit)
#
link.se=suppressWarnings(sqrt(chat*diag(design%*%model$vcv[indices,indices]%*%t(design))))
link.se[is.na(link.se)]=0
links=rep(link,length(real))
ind=unique(c(grep("mlogit",links,ignore.case=TRUE),which(links%in%c("sin","Sin","LogLog","loglog","CLogLog","cloglog"))))
linkse=suppressWarnings(sqrt(diag(vcv.real)[ind])/(real[ind]*(1-real[ind])))
linkse[is.na(linkse)]=0
linkse[is.infinite(linkse)]=0
link.se[ind]=linkse
link.values=design%*%beta
link.values[ind]=suppressWarnings(log(real[ind]/(1-real[ind])))
link.values[ind][abs(real[ind]-1)<1e-7]=100
link.values[ind][abs(real[ind]-0)<1e-7]=-100
links[ind]="logit"
real.lcl=convert.link.to.real(link.values-1.96*link.se,links=links)
real.ucl=convert.link.to.real(link.values+1.96*link.se,links=links)
#
# Set v-c values of fixed parameters to 0
#
vcv.real[fixedparms,]=0
vcv.real[,fixedparms]=0
diag(vcv.real)[diag(vcv.real)<0]=0
se.real=sqrt(diag(vcv.real))
#se.real[is.na(se.real)]=0
fixed=rep("",dim(design)[1])
fixed[fixedparms]="Fixed"
#fixed[boundaryparms]="Boundary"
if(vcv)
   return(list(real=cbind(df,data.frame(estimate=real,se=se.real,lcl=real.lcl,ucl=real.ucl,fixed=fixed)),vcv.real=vcv.real))
else
   return(cbind(df,data.frame(estimate=real,se=se.real,lcl=real.lcl,ucl=real.ucl,fixed=fixed)))
}

