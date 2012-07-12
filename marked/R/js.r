#' Fitting function for Jolly-Seber model using Schwarz-Arnason POPAN
#' formulation
#' 
#' A function for computing MLEs for a specified Jolly-Seber open population
#' capture-recapture model for processed dataframe \code{x} with user specified
#' formulas in \code{parameters} that create list of design matrices
#' \code{dml}. This function can be called directly but is most easily called
#' from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{js} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' Be cautious with this function at present.  It does not include many checks
#' to make sure values like fixed values will remain in the specified range of
#' the data.  Normally this would not be a big problem but because
#' \code{\link{js.lnl}} calls an external FORTRAN subroutine via
#' \code{\link{cjs.lnl}}, if it gets a subscirpt out of bounds, it will cause R
#' to terminate.  So make sure to save your workspace frequently if you use
#' this function in its current implementation.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, Phi.dm,p.dm,Phi.fixed,p.fixed, and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to js when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for Phi and p to speed up execution
#' @param Phi initial value of Phi; used to set intercept parameter
#' @param p initial value of p; used to set intercept parameter
#' @param initial initial values for parameters if desired; if named vector
#' from previous run it will match to columns with same name
#' @param method method to use for optimization; see \code{optimx}
#' @param hessian if TRUE will compute and return the hessian
#' @param debug if TRUE will print out information for each iteration
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param control control string for optimization functions
#' @param scale vector of scale values for parameters
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optimx} or \code{\link{js.lnl}}
#' @return The resulting value of the function is a list with the class of
#' crm,js such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optimx}; if 0 \code{optimx} thinks it
#' converged} \item{count}{\code{optimx} results of number of function
#' evaluations} \item{reals}{dataframe of data and real Phi and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Schwarz, C. J., and A. N. Arnason. 1996. A general methodology
#' for the analysis of capture-recapture experiments in open populations.
#' Biometrics 52:860-873.
js=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,Phi=NULL,p=NULL,initial=NULL,method="BFGS",
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,...)
{
#
#  Setup values from arguments
#    Phi.dm          - Phi design matrix created by call to create.dm
#    p.dm            - p design matrix created by call to create.dm
#    pent.dm         - pent design matrix created by call to create.dm
#    N.dm            - N design matrix created by call to create.dm
#    Phi.dmdf        - Phi design dataframe created by call to create.dmdf
#    p.dmdf          - p design dataframe created by call to create.dmdf
#    pent.dmdf       - pent design dataframe created by call to create.dmdf
#    N.dmdf          - N design dataframe created by call to create.dmdf
#    Phi.fixed       - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix phi(i,j)=f 
#    p.fixed         - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
#    pent.fixed      - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix pent(i,j)=f 
   Phi.dmdf=ddl$Phi
   p.dmdf=ddl$p
   nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal
   time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
   if(!is.null(Phi.dmdf$time.interval))
	   time.intervals=matrix(Phi.dmdf$time.interval,nrow(x$data),ncol=nocc-1,byrow=T)
# Compute nobstot - total unique critters
   if(nrow(dml$N)==1) 
      nobstot=sum(x$data$freq)    
   else
     nobstot=tapply(x$data$freq,x$data$group,sum)   
#  If no fixed real parameters are specified, assign dummy unused ones with negative indices and 0 value
   if(is.null(parameters$Phi.fixed))
	   parameters$Phi.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)
   if(is.null(parameters$p.fixed))
	   parameters$p.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
   if(is.null(parameters$pent.fixed))
	   parameters$pent.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
#  Store data from x into x
   x=x$data
#  set default frequencies if not used
   if(is.null(x$freq))
	   freq=NULL
   else
	   freq=x$freq
 #  get first and last vectors, loc and chmat
   ch=x$ch
   imat=process.ch(ch,freq,all=TRUE)
#  Use specified initial values or create if null
   if(is.null(initial))
   {
	   if(is.null(Phi)|is.null(p)) 
		   par=cjs.initial(dml,imat)
	   else  
		   par=c(log(Phi/(1-Phi)),rep(0,ncol(dml$Phi)-1),log(p/(1-p)),rep(0,ncol(dml$p)-1))
	   par=c(par,0,rep(0,ncol(dml$pent)-1),1,rep(0,ncol(dml$N)-1))
   } 
   else
   {
	   if(is.null(names(initial)))
	   {
		   if(length(initial)!=(ncol(dml$Phi)+ncol(dml$p)+ncol(dml$pent)+ncol(dml$N)))
			   stop("Length of initial vector does not match number of parameters.")
		   else
			   par=initial
	   }
	   else
	   {
		   beta.names=c(paste("Phi:",colnames(dml$Phi),sep="") ,paste("p:",colnames(dml$p),sep=""),paste("pent:",colnames(dml$pent),sep=""),paste("N:",colnames(dml$N),sep=""))
		   par=rep(0,length(beta.names))
		   par[beta.names%in%names(initial)]=initial[which(names(initial)%in%beta.names)]
	   }
   }
#  Create list of model data for optimization; if passed as an argument create model_data.save 
#  and use model_data (accumulated values); otherwise create model_data, save it and accumulate it
#  if requested.
   if(!is.null(model_data)) 
   {
	   model_data.save=list(Phi.dm=dml$Phi,p.dm=dml$p,pent.dm=dml$pent,N.dm=dml$N,imat=imat,Phi.fixed=parameters$Phi.fixed,
			   p.fixed=parameters$p.fixed,pent.fixed=parameters$pent.fixed,time.intervals=time.intervals)
   }else
   {
	   model_data=list(Phi.dm=dml$Phi,p.dm=dml$p,pent.dm=dml$pent,N.dm=dml$N,imat=imat,Phi.fixed=parameters$Phi.fixed,
			   p.fixed=parameters$p.fixed,pent.fixed=parameters$pent.fixed,time.intervals=time.intervals)
#       If data are to be accumulated based on ch and design matrices do so here;
	   if(accumulate)
	   {
		   cat("\n Accumulating capture frequencies based on design. This can take awhile.\n")
		   flush.console()
		   model_data.save=model_data   
		   model_data=js.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	   }else
		   model_data.save=NULL
   }
   
   #  Scale the design matrices and parameters with either input scale or computed scale
   if(is.null(scale))
   {
	   scale.phi=apply(model_data$Phi.dm,2,function(x) mean(x[x!=0]))
	   scale.p=apply(model_data$p.dm,2,function(x) mean(x[x!=0]))
	   scale.pent=apply(model_data$pent.dm,2,function(x) mean(x[x!=0]))
	   scale.N=apply(model_data$N.dm,2,function(x) mean(x[x!=0]))
   } else
   {
	   if(all(scale==1))
	   {
		   scale.phi=rep(1,ncol(model_data$Phi.dm))
		   scale.p=rep(1,ncol(model_data$p.dm))
		   scale.pent=rep(1,ncol(model_data$pent.dm))
		   scale.N=rep(1,ncol(model_data$N.dm))
	   }else
	   {
		   len=ncol(model_data$Phi.dm)
		   scale.phi=scale[1:len]
		   scale.p=scale[(len+1):(len+ncol(model_data$p.dm))]
		   len=len+ncol(model_data$p.dm)
		   scale.pent=scale[(len+1):(len+ncol(model_data$pent.dm))]
		   len=len+ncol(model_data$pent.dm)
		   scale.N=scale[(len+1):length(scale)]
	   }
   }
   model_data$Phi.dm=Matrix:::t(Matrix:::t(model_data$Phi.dm)/scale.phi)
   model_data$p.dm=Matrix:::t(Matrix:::t(model_data$p.dm)/scale.p)
   model_data$pent.dm=Matrix:::t(Matrix:::t(model_data$pent.dm)/scale.pent)
   model_data$N.dm=Matrix:::t(Matrix:::t(model_data$N.dm)/scale.N)
   par=par*c(scale.phi,scale.p,scale.pent,scale.N)
#  call optim to find mles with js.lnl which gives -log-likelihood
   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
   cat("\n Starting optimization",ncol(model_data$Phi.dm)+ncol(model_data$p.dm)+ncol(model_data$pent.dm)+ncol(model_data$N.dm)," parameters\n")
   convergence=1
   i=0
   while (convergence!=0 & i <= refit)
   {
	   if(i>0)
	   {
		   cat("\n Re-starting optimization using Nelder-Mead for ",ncol(model_data$Phi.dm)+ncol(model_data$p.dm)+ncol(model_data$pent.dm)+ncol(model_data$N.dm)," parameters\n")   
		   method="Nelder-Mead"
		   itnmax=2*itnmax
	   }
	   mod=suppressPackageStartupMessages(optimx(par,js.lnl,model_data=model_data,method=method,hessian=hessian,
					   debug=debug,control=control,itnmax=itnmax,nobstot=nobstot,...))
	   par=mod$par$par
	   convergence=mod$conv$conv
	   i=i+1
	   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
   }
   js.beta=mod$par$par/c(scale.phi,scale.p,scale.pent,scale.N)
   names(js.beta)=c(paste("Phi:",colnames(model_data$Phi.dm),sep="") ,paste("p:",colnames(model_data$p.dm),sep=""),paste("pent:",colnames(model_data$pent.dm),sep=""),paste("N:",colnames(model_data$N.dm),sep=""))
   if(is.null(model_data.save))
   {
	  if(is.null(x$group))
		  ui=tapply(model_data$imat$freq,list(model_data$imat$first),sum)
	  else
		  ui=tapply(model_data$imat$freq,list(model_data$imat$first,x$group),sum)	  
   } else
   {
      if(is.null(x$group))
	      ui=tapply(model_data.save$imat$freq,list(model_data.save$imat$first),sum)
      else
	      ui=tapply(model_data.save$imat$freq,list(model_data.save$imat$first,x$group),sum)
   }
   lnl=mod$fvalues$fvalues+sum(lfactorial(ui))
   res=list(beta=js.beta,neg2lnl=2*lnl,AIC=2*lnl+2*length(js.beta),
		   convergence=mod$conv,count=mod$itns,
		   scale=list(phi=scale.phi,p=scale.p,pent=scale.pent,N=scale.N),
		   model_data=model_data)
#  Restore complete non-accumulated model_data with unscaled design matrices and compte reals
   if(!is.null(model_data.save)) model_data=model_data.save
#  compute reals
   nphi=ncol(model_data$Phi.dm)
   np=ncol(model_data$p.dm)
   npent=ncol(model_data$pent.dm)
   nN=ncol(model_data$N.dm)
   beta.phi=js.beta[1:nphi]
   beta.p=js.beta[(nphi+1):(nphi+np)]
   beta.pent=js.beta[(nphi+np+1):(nphi+np+npent)]
   beta.N=js.beta[(nphi+np+npent+1):(nphi+np+npent+nN)]
# create Phi and p beta matrices excluding first occasion on p
   Phis=plogis(as.vector(model_data$Phi.dm%*%beta.phi))
   ps=plogis(as.vector(model_data$p.dm%*%beta.p))
   pents=cbind(rep(1,nrow(model_data$pent.dm)),exp(as.vector(model_data$pent.dm%*%beta.pent)))
   pents=pents/apply(pents,1,sum)
   Ns=exp(as.vector(model_data$N.dm%*%beta.N))
#  
   reals=p.dmdf
   names(reals)[names(reals)=="time"]="p.time"
   names(reals)[names(reals)=="age"]="p.age"
   Phirecs=Phi.dmdf[,!colnames(Phi.dmdf)%in%colnames(reals),drop=FALSE]
   Phirecs[(rep(1:nrow(x),each=nocc)-1)*(nocc-1)+c(1:(nocc-1),nocc-1),]
   Phirecs[seq(nocc,nrow(x)*nocc,by=nocc),]=NA
   names(reals)[names(reals)=="time"]="Phi.time"
   names(reals)[names(reals)=="age"]="Phi.age"
   res$reals=reals
# replace with call to js.hessian#  If requested compute hessian and var-cov matrix 
   if(hessian) 
   {
	   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
	   cat("\n Computing hessian\n")
	   res$vcv=js.hessian(res, nobstot)
   }   
   class(res)=c("crm","js")
   return(res)
}


