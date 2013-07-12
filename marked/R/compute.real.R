#' Compute estimates of real parameters
#' 
#' Computes real estimates and their var-cov for a particular subset of 
#' parameters.
#' 
#' @param model model object
#' @param parameter name of real parameter to be computed (eg "Phi" or "p")
#' @param ddl list of design data 
#' @param unique TRUE if only unique values should be returned
#' @param vcv logical; if TRUE, computes and returns v-c matrix of real estimates
#' @param se logical; if TRUE, computes std errors and conf itervals of real estimates
#' @param chat over-dispersion value
#' @param subset logical expression using fields in real dataframe
#' @param select character vector of field names in real that you want to include 
#' @param showDesign if TRUE, show design matrix instead of data 
#' @export
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
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~sex+time),p=list(formula=~1)),hessian=TRUE)
#' xx=compute.real(mod.Phisex.pdot,"Phi",unique=TRUE,vcv=TRUE)
#' @keywords utility
compute.real <-function(model,parameter,ddl=NULL,unique=TRUE,vcv=FALSE,se=FALSE,chat=1,subset,select,showDesign=FALSE)
{
  if(is.null(ddl))return(model$results$reals)
  mcmc=ifelse(class(model)[2]=="mcmc",TRUE,FALSE)
# Setup data and dm
  data=ddl[[parameter]]
  design=model$results$model_data[[paste(parameter,"dm",sep=".")]]
  # remove unneeded reals
  if(!is.null(data$Time)&model$model.parameter[[parameter]]$type=="Triang")
  {
	  indices=which(data$Time>=data$Cohort)
	  data=data[indices,,drop=FALSE]
	  design=design[indices,,drop=FALSE]
  }
  if(parameter!="N" &!is.null(data$freq))
  {
	  indices=which(data$freq!=0)
	  data=data[indices,,drop=FALSE]
	  design=design[indices,,drop=FALSE]
  }
  design=as.matrix(design)
  # Set up links
  link=model$model.parameter[[parameter]]$link
  if(link=="mlogit")
	  if(parameter=="pent")
		  link=paste("mlogit",data$id,sep="")
  if(length(link)>1)link=link[indices]
  # set up data to be displayed with estimates
  df=data
  if(missing(select))
  {
	  if(unique)
	      {df=df[,all.vars(model$model.parameters[[parameter]]$formula),drop=FALSE]}
  } else
      df=df[,select,drop=FALSE]	  
  # Check to make sure beta and DM match
  results=model$results
  beta=results$beta[[parameter]]
  if(mcmc)
  {
	  if(ncol(design)!=ncol(model$results$beta.mcmc[[parameter]]))
   	     stop("Mismatch between number of design columns and length of beta")
  } else
  {
	  if(ncol(design)!=length(beta))
         stop("Mismatch between number of design columns and length of beta")
  }
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
   if(mcmc)
	   real=convert.link.to.real(t(design%*%t(results$beta.mcmc[[parameter]])),links=link,fixed=fixedvalues)
   else
       real=convert.link.to.real(design%*%beta,links=link,fixed=fixedvalues)
#
#  Set fixed real parameters to their fixed values
#
   real[fixedparms]=fixedvalues[fixedparms]
#
#  Summarize reals for mcmc models
#
if(mcmc)
{
    fitted <- mcmc(real)
    summ <- summary(fitted)
    hpd <- HPDinterval(fitted)
    reals <- data.frame(
		  mode = apply(real, 2, mcmc_mode), 
  	      mean=apply(real, 2, mean), 
          sd = apply(real, 2, sd),
		  CI.lower=hpd[,1], CI.upper=hpd[,2])
} else
#  If no vcv requested only create reals dataframe
   if(!vcv&!se)
   {
	   reals=data.frame(estimate=real)
	   rownames(reals)=NULL
   } else
   {
	   if(is.null(results$beta.vcv))stop("vcv matrix not available in model")
#      Compute vc matrix for real parameters which needs to be modified for
#      any mlogit parameters
	   indices=grep(paste(parameter,"\\.",sep=""),colnames(results$beta.vcv))
	   if(any(substr(link,1,6)!="mlogit"))
	   {
		   deriv.real=as.matrix(deriv_inverse.link(real,design,link))
		   vcv.real=deriv.real%*%results$beta.vcv[indices,indices,drop=FALSE]%*%t(deriv.real)
	   } else
#      If vcv=TRUE, compute v-c matrix and std errors of real estimates
#      To handle any mlogit parameters compute pseudo-real estimates using log in place of mlogit
	   {
		   pseudo.real=as.vector(convert.link.to.real(design%*%beta,links="log"))
		   pseudo.real[fixedparms]=fixedvalues[fixedparms]
		   pseudo.real[fixedparms]=exp(pseudo.real[fixedparms])
#          Compute first derivatives of pseudo-real (for any mlogit parameter)
#          estimates with respect to beta parameters
		   deriv.pseudo=deriv_inverse.link(pseudo.real,design,"log")
		   deriv.pseudo[fixedparms,]=0
		   vcv.pseudo=deriv.pseudo%*%results$beta.vcv[indices,indices,drop=FALSE]%*%t(deriv.pseudo)
#          Apply chain rule to get variance of real parameters which has mlogits
#          expressed as zi/(1+z1+...zk) where k is number of mlogit components-1 and
#          non-mlogits are expressed as zi.
#          bottom is either 1 for non-mlogits and the sum for mlogits
#          pbottom is partial with respect to zi
		   links=link
		   mlogits=outer(links,links,function(x,y)as.numeric(x==y))*as.numeric(substr(links,1,6)=="mlogit")
		   pbottom=matrix(0,nrow=dim(vcv.pseudo)[1],ncol=dim(vcv.pseudo)[1]) + mlogits
		   bottom=diag(nrow=nrow(vcv.pseudo))*(1-as.numeric(substr(links,1,6)=="mlogit"))+
				   mlogits + pbottom*apply(pbottom*pseudo.real,2,sum)
		   deriv.pseudo=(diag(nrow=dim(vcv.pseudo)[1])*bottom-pseudo.real*pbottom)/bottom^2
		   deriv.pseudo[is.nan(deriv.pseudo)]=0
		   vcv.real=deriv.pseudo%*%vcv.pseudo%*%t(deriv.pseudo)
	   }
	   vcv.real=chat*vcv.real
#      Compute conf interval taking into account use of logit transform for mlogit
#      and any 0-1 link (loglog,cloglog,sin,logit)
	   link.se=suppressWarnings(sqrt(chat*diag(design%*%results$beta.vcv[indices,indices]%*%t(design))))
	   link.se[is.na(link.se)]=0
	   if(length(link)==1)
		   links=rep(link,length(real))
	   else
		   links=link
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
#      Set v-c values of fixed parameters to 0
	   vcv.real[fixedparms,]=0
	   vcv.real[,fixedparms]=0
	   diag(vcv.real)[diag(vcv.real)<0]=0
	   se.real=sqrt(diag(vcv.real))
#se.real[is.na(se.real)]=0
	   fixed=rep("",dim(design)[1])
	   fixed[fixedparms]="Fixed"
	   reals=data.frame(estimate=real,se=se.real,lcl=real.lcl,ucl=real.ucl)    
#fixed[boundaryparms]="Boundary"   
    }   
#  Setup subset or get unique records 
    if(!missing(subset))
	{
		if(subset!="")
		{
			if(substitute(subset)!="substitute(subset)")
				sub=substitute(subset)
			else
				sub=subset
			indices = eval(sub, data)
		}
		else
		{
			if(unique)
				indices=which(!duplicated(apply(design,1,paste,collapse="")))
			else
				indices = TRUE
		}
	}
	df=df[indices,,drop=FALSE]	  
	design=design[indices,,drop=FALSE]
	if(ncol(df)==0)df$Intercept=1
    reals=reals[indices,,drop=FALSE]
	if(showDesign)
		reals=cbind(design,reals)
	else
		reals=cbind(df,reals)
	reals=reals[do.call(order, reals),]
	rownames(reals)=NULL
	if(vcv)
		return(list(real=reals,vcv=vcv.real[indices,indices]))
    else
        return(reals)
}
# Temp function for HMM to be merged into compute.real at some point
# computes real estimates using inverse of link from design data (ddl) and model for a particular parameter type (parname) or
# returns the number of columns in the design matrix (compute=FALSE); handles fixed parameters assigned by non-NA value in field named 
# fix in the ddl dataframe.
reals=function(parname,ddl,parameters,parlist=NULL,compute=TRUE)
{
	# create design matrix (dm) for parameter parname
	dm=model.matrix(parameters[[parname]]$formula,ddl)
	# if some reals are fixed, assign 0 to rows of dm and then
	# remove any columns (parameters) that are all 0.
	if(!is.null(ddl$fix))
	{
		dm[!is.na(ddl$fix),]=0
		dm=dm[,apply(dm,2,function(x) any(x!=0)),drop=FALSE]
	}
	# if not computing reals, return the number of colmns in dm
	if(!compute)return(ncol(dm))
	# Currently for log or logit link, return the inverse values
	if(parameters[[parname]]$link=="log")
		values=exp(dm%*%parlist[[parname]])
	else if(parameters[[parname]]$link=="logit")
		values=plogis(dm%*%parlist[[parname]])
	# if some reals are fixed, set reals to their fixed values 
	if(!is.null(ddl$fix))
		values[!is.na(ddl$fix)]=ddl$fix[!is.na(ddl$fix)]
	# return vector of reals
	return(values)
}
