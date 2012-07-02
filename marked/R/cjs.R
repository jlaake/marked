#' Fitting function for CJS models
#' 
#' A function for computing MLEs for a specified Cormack-Jolly-Seber open
#' population capture-recapture model for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{cjs} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' Be cautious with this function at present.  It does not include many checks
#' to make sure values like fixed values will remain in the specified range of
#' the data.  Normally this would not be a big problem but because
#' \code{\link{cjs.lnl}} calls an external FORTRAN subroutine, if it gets a
#' subscript out of bounds, it will cause R to terminate.  So make sure to save
#' your workspace frequently if you use this function in its current
#' implementation.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, Phi.dm,p.dm,Phi.fixed,p.fixed, and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for Phi and p to speed up execution
#' @param Phi initial value of Phi; used to set intercept parameter
#' @param p initial value of p; used to set intercept parameter
#' @param initial initial values for parameters if desired; if named vector
#' from previous run it will match to columns with same name
#' @param method method to use for optimization; see \code{optim}
#' @param hessian if TRUE will compute and return the hessian
#' @param debug if TRUE will print out information for each iteration
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param control control string for optimization functions
#' @param scale vector of scale values for parameters
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optim} or \code{\link{cjs.lnl}}
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optim}; if 0 \code{optim} thinks it
#' converged} \item{count}{\code{optim} results of number of function
#' evaluations} \item{reals}{dataframe of data and real Phi and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laake <jeff.laake@@noaa.gov>
cjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,Phi=NULL,p=NULL,initial=NULL,method,
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale, ...)
{
#
#  Setup values from arguments
#    Phi.dm          - Phi design matrix created by call to create.dm
#    p.dm            - p design matrix created by call to create.dm
#    Phi.dmdf        - Phi design dataframe created by call to create.dmdf
#    p.dmdf          - p design dataframe created by call to create.dmdf
#    Phi.fixed       - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix phi(i,j)=f 
#    p.fixed         - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
   Phi.dmdf=ddl$Phi
   p.dmdf=ddl$p
   nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use default of 1 if none are in Phi.dmdf
   time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
   if(!is.null(Phi.dmdf$time.interval))
	   time.intervals=matrix(Phi.dmdf$time.interval,nrow(x$data),ncol=nocc-1,byrow=T)
#  If no fixed real parameters are specified, assign dummy unused ones with negative indices and 0 value
   if(is.null(parameters$Phi$fixed))
	   parameters$Phi$fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)
   if(is.null(parameters$p$fixed))
	   parameters$p$fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
#  Store data from x into x
   x=x$data
#  set default frequencies if not used
   if(is.null(x$freq))
      freq=NULL
   else
      freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
   ch=x$ch
   imat=process.ch(ch,freq)
#  Use specified initial values or create if null
   if(is.null(initial))
   {
	   if(is.null(Phi)|is.null(p))
		   par=cjs.initial(dml,imat)
	   else  
	       par=c(log(Phi/(1-Phi)),rep(0,ncol(dml$Phi)-1),log(p/(1-p)),rep(0,ncol(dml$p)-1))
   }
   else
   {
	   if(is.null(names(initial)))
	   {
		   if(length(initial)!=(ncol(dml$Phi)+ncol(dml$p)))
			   stop("Length of initial vector does not match number of parameters.")
		   else
			   par=initial
	   }
	   else
	   {
		   beta.names=c(paste("Phi:",colnames(dml$Phi),sep="") ,paste("p:",colnames(dml$p),sep=""))
		   par=rep(0,length(beta.names))
		   par[beta.names%in%names(initial)]=initial[which(names(initial)%in%beta.names)]
	   }
   }
#  Create list of model data for optimization; if passed as an argument create model_data.save 
#  and use model_data (accumulated values); otherwise create model_data, save it and accumulate it
#  if requested.
    if(!is.null(model_data)) 
	{
		model_data.save=list(Phi.dm=dml$Phi,p.dm=dml$p,imat=imat,Phi.fixed=parameters$Phi$fixed,
				p.fixed=parameters$p$fixed,time.intervals=time.intervals)
	}else
	{
		model_data=list(Phi.dm=dml$Phi,p.dm=dml$p,imat=imat,Phi.fixed=parameters$Phi$fixed,
				p.fixed=parameters$p$fixed,time.intervals=time.intervals)
#       If data are to be accumulated based on ch and design matrices do so here;
		if(accumulate)
		{
			cat("\n Accumulating capture frequencies based on design. This can take awhile.\n")
			flush.console()
			model_data.save=model_data   
			model_data=cjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
		}else
			model_data.save=NULL
	}
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#
#  Scale the design matrices and parameters with either input scale or computed scale
   if(is.null(scale))
   {
      scale.phi=apply(model_data$Phi.dm,2,function(x) mean(x[x!=0]))
	  scale.p=apply(model_data$p.dm,2,function(x) mean(x[x!=0]))
   } else
   {
	   if(all(scale==1))
	   {
		   scale.phi=rep(1,ncol(model_data$Phi.dm))
		   scale.p=rep(1,ncol(model_data$p.dm))
	   }else
	   {
		   scale.phi=scale[1:ncol(model_data$Phi.dm)]
		   scale.p=scale[(ncol(model_data$Phi.dm)+1):length(scale)]
	   }
   }
   model_data$Phi.dm=Matrix:::t(Matrix:::t(model_data$Phi.dm)/scale.phi)
   model_data$p.dm=Matrix:::t(Matrix:::t(model_data$p.dm)/scale.p)
   par=par*c(scale.phi,scale.p)
#  Call optimx to find mles with cjs.lnl which gives -2 * log-likelihood
   cat("\n Starting optimization for ",ncol(model_data$Phi.dm)+ncol(model_data$p.dm)," parameters\n")
   flush.console()
   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
   convergence=1
   i=0
   while (convergence!=0 & i <= refit)
   {
       if(i>0)
	   {
		   cat("\n Re-starting optimization using Nelder-Mead for ",ncol(model_data$Phi.dm)+ncol(model_data$p.dm)," parameters\n")   
		   method="Nelder-Mead"
           itnmax=2*itnmax
	   }
	   mod=suppressPackageStartupMessages(optimx(par,cjs.lnl,model_data=model_data,Phi.links=NULL,p.links=NULL,method=method,hessian=FALSE,
			         debug=debug,control=control,itnmax=itnmax,...))
	  par=mod$par$par
	  convergence=mod$conv$conv
	  i=i+1
	  assign(".markedfunc_eval", 0, envir = .GlobalEnv)
  }
#  Rescale parameter vector, restore model_data and call cjs.lnl to compute all values and not just -2lnl
   cjs.beta=par/c(scale.phi,scale.p)
   names(cjs.beta)=c(paste("Phi:",colnames(model_data$Phi.dm),sep="") ,paste("p:",colnames(model_data$p.dm),sep=""))
#  Create results list 
   lnl=mod$fvalues$fvalues
   res=list(beta=cjs.beta,neg2lnl=lnl,AIC=lnl+2*length(cjs.beta),convergence=mod$conv,count=mod$itns,mod=mod,scale=list(phi=scale.phi,p=scale.p),model_data=model_data)
#  Restore complete non-accumulated model_data with unscaled design matrices and call cjs to compute all values including reals
   if(!is.null(model_data.save)) model_data=model_data.save
   allval=cjs.lnl(cjs.beta,model_data.save,Phi.links=NULL,p.links=NULL,all=TRUE)
#  Create dataframe of real parameter estimates
   reals=Phi.dmdf
   names(reals)[names(reals)=="time"]="Phi.time"
   names(reals)[names(reals)=="age"]="Phi.age"
   reals=cbind(reals,p.dmdf[,!colnames(p.dmdf)%in%colnames(reals)])
   names(reals)[names(reals)=="time"]="p.time"
   names(reals)[names(reals)=="age"]="p.age"
   xx=matrix(allval[[2]],nrow=nrow(x),ncol=nocc-1)
   reals$Phi=as.vector(t(xx))
   xx=matrix(allval[[3]],nrow=nrow(x),ncol=nocc-1)
   reals$p=as.vector(t(xx))
   res$reals=reals[reals$Time>=reals$Cohort,]
#  If requested compute hessian and var-cov matrix 
   if(hessian) 
   {
	   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
	   cat("\n Computing hessian\n")
	   res$vcv=cjs.hessian(res, Phi.links=NULL, p.links=NULL)
   }   
#  Assign S3 class and return
   class(res)=c("crm","cjs")
   return(res)
}


