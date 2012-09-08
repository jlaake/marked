#' Compute estimates of real parameters
#' 
#' Computes real estimates and their var-cov for a particular subset of 
#' parameters.
#' 
#' @param model model object
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
#' @examples
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~sex+time),p=list(formula=~1)),hessian=T)
#' xx=compute.real(mod.Phisex.pdot,"Phi",unique=T,vcv=TRUE)
#' @keywords utility
compute.real <-function(model,parameter,unique=TRUE,vcv=FALSE,se=FALSE,chat=1,subset,select)
{
  links=c(Phi="logit",p="logit",pent="mlogit",N="log")
  link=links[parameter]
  if(!missing(subset))
  {
	  sub = substitute(subset)
	  indices = eval(sub, model$results$reals)
  }
  else
	  indices = TRUE
  if(missing(select))
  {
	  design=model$design.parameters[[parameter]]
	  realnames=names(model$results$reals)
	  parspecific=grep(paste(parameter,".",sep=""),realnames)
	  if(length(parspecific))
	     realnames=realnames[parspecific]
	  else
		 realnames=NULL
	  select=c("cohort",realnames,design$static,design$time.varying)
  }
  df=model$results$reals[indices,select,drop=FALSE]
 # get dm
	design=model$results$model_data[[paste(parameter,"dm",sep=".")]][model$results$reals$rec,,drop=FALSE]
	design=design[indices,,drop=FALSE]
# if unique, find unique dm values
  if(unique)
  {
	  indices=which(!duplicated(apply(design,1,paste,collapse="")))
	  df=df[indices,,drop=FALSE]
	  design=design[indices,,drop=FALSE]
  }
  design=as.matrix(design)
#
# Check to make sure beta and DM match
#
beta=model$results$beta[grep(paste(parameter,":",sep=""),names(model$results$beta))]
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
if(!vcv&!se)return(cbind(df,data.frame(estimate=real)))
if(is.null(model$results$beta.vcv))stop("vcv matrix not available in model")
#
# Compute vc matrix for real parameters which needs to be modified for
# any mlogit parameters
#
indices=grep(paste(parameter,":",sep=""),colnames(model$results$beta.vcv))
if(link!="mlogit")
{
	deriv.real=as.matrix(deriv_inverse.link(real,design,link))
	vcv.real=deriv.real%*%model$results$beta.vcv[indices,indices,drop=FALSE]%*%t(deriv.real)
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
    vcv.pseudo=deriv.pseudo%*%model$results$beta.vcv[indices,indices,drop=FALSE]%*%t(deriv.pseudo)
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
link.se=suppressWarnings(sqrt(chat*diag(design%*%model$results$beta.vcv[indices,indices]%*%t(design))))
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
   return(list(real=cbind(df,data.frame(estimate=real,se=se.real,lcl=real.lcl,ucl=real.ucl)),vcv.real=vcv.real))
else
   return(cbind(df,data.frame(estimate=real,se=se.real,lcl=real.lcl,ucl=real.ucl)))
}

