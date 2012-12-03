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
#' @param initial list of initial values for parameters if desired; if each is a named vector
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
#' @param use.admb if TRUE creates data file for admbcjs.tpl and runs admb optimizer
#' @param re if TRUE creates random effect model admbcjsre.tpl and runs admb optimizer
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optim} or \code{\link{cjs.lnl}}
#' @import R2admb
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
#' @references Pledger, S., K. H. Pollock, et al. (2003). Open
#' capture-recapture models with heterogeneity: I. Cormack-Jolly-Seber model.
#' Biometrics 59(4):786-794.
cjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
			use.admb=FALSE,re=FALSE,compile=FALSE,extra.args="",...)
{
   if(re)accumulate=FALSE
   nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use default of 1 if none are in Phi.dmdf
   time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
   if(!is.null(ddl$Phi$time.interval))
	   time.intervals=matrix(ddl$Phi$time.interval,nrow(x$data),ncol=nocc-1,byrow=TRUE)
#  If no fixed real parameters are specified, assign dummy unused ones with negative indices and 0 value
   if(is.null(parameters$Phi$fixed))
	   parameters$Phi$fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)
   if(is.null(parameters$p$fixed))
	   parameters$p$fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
#  Store data from x into x
   x=x$data
#  set default frequencies if not used
   freq=NULL
   if(!is.null(x$freq))freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
   ch=x$ch
   imat=process.ch(ch,freq,all=FALSE)
#  Use specified initial values or create if null
   if(is.null(initial))
	   par=cjs.initial(dml,imat)
   else
       par=set.initial(names(dml),dml,initial)
#  Create list of model data for optimization
	model_data=list(Phi.dm=dml$Phi,p.dm=dml$p,imat=imat,Phi.fixed=parameters$Phi$fixed,
			p.fixed=parameters$p$fixed,time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
	if(accumulate)
	{
		cat("Accumulating capture histories based on design. This can take awhile.\n")
		flush.console()
		model_data.save=model_data   
		model_data=cjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#  Scale the design matrices and parameters with either input scale or computed scale
   if(use.admb)scale=1
   scale=set.scale(names(dml),model_data,scale)
   model_data=scale.dm(model_data,scale)
   if(re)use.admb=TRUE
   if(!use.admb)
   {
	   par=scale.par(par,scale)
	   #  Call optimx to find mles with cjs.lnl which gives -log-likelihood
	   cat("Starting optimization for ",length(par)," parameters\n")
	   flush.console()
	   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
	   if("SANN"%in%method)
	   {
		   mod=optim(par,cjs.lnl,model_data=model_data,method="SANN",hessian=FALSE,
				   debug=debug,control=control,...)
		   par= mod$par
		   convergence=mod$convergence
		   lnl=mod$value
		   counts=mod$counts
	   }else
	   {
		   mod=suppressPackageStartupMessages(optimx(par,cjs.lnl,model_data=model_data,method=method,hessian=FALSE,
						   debug=debug,control=control,itnmax=itnmax,...))
		   objfct=unlist(mod$fvalues)
		   bestmin=which.min(objfct)
		   par= mod$par[[bestmin]]
		   convergence=mod$conv[[bestmin]]
		   counts=mod$itns[[length(mod$itns)]]
		   lnl=mod$fvalues[[bestmin]]
	   }
	   #  Rescale parameter vector 
	   cjs.beta=unscale.par(par,scale)
       # Create results list 
	   res=list(beta=cjs.beta,neg2lnl=2*lnl,AIC=2*lnl+2*sum(sapply(cjs.beta,length)),
			   convergence=convergence,count=counts,optim.details=mod,
			   scale=scale,model_data=model_data,
			   options=list(accumulate=accumulate,initial=initial,method=method,
					   chunk_size=chunk_size,itnmax=itnmax,control=control))
       # Compute hessian if requested
	   if(hessian) 
	   {
		   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
		   cat("Computing hessian\n")
		   res$beta.vcv=cjs.hessian(res)
	   } 
	   assign(".markedfunc_eval", 0, envir = .GlobalEnv)	   
   } else
   {
       # see if admb can be found; this is not a complete test but should catch the novice user who has
	   # not setup admb at all
	   if(Sys.which("tpl2cpp.exe")=="")stop("admb not found; setup links to admb and c++ compiler with environment variables or put in path") 
       # cleanup any leftover admbcjs files
	   sdir=system.file(package="marked")
	   # if admbcjs.tpl is not available copy from the package directory
	   if(!re)
	      tpl="admbcjs"
       else
	   {
		   tpl="admbcjsre"
		   if(!hessian)message("ignoring hessian setting; set to TRUE")
		   hessian=TRUE
	   }
	   clean_admb(tpl)
	   if(!file.exists(paste(tpl,".tpl",sep="")))
		   file.copy(file.path(sdir,paste(tpl,".tpl",sep="")),file.path(getwd(),paste(tpl,".tpl",sep="")),overwrite=TRUE)
	   if(!file.exists(paste(tpl,".exe",sep=""))|compile)
		   compile_admb(tpl,re=re)
	   # create admbcjs.dat file to create its contents 
	   con=file(paste(tpl,".dat",sep=""),open="wt")
	   # Number of observations
	   n=length(model_data$imat$freq)
	   write(n,con,append=FALSE)
	   # Number of occasions
	   nocc=model_data$imat$nocc
	   write(nocc,con,append=TRUE)
	   # capture history matrix
	   write(t(model_data$imat$chmat),con,ncolumns=nocc,append=TRUE)
	   # first occasions seen 
	   write(model_data$imat$first,con,ncolumns=n,append=TRUE)
	   # last occasions seen 
	   write(model_data$imat$last,con,ncolumns=n,append=TRUE)
	   # frequency of capture history 
	   if(!re)
	   {
	      write(model_data$imat$freq,con,ncolumns=n,append=TRUE)
	   } else
	   {
		   if(any(model_data$imat$freq!=1))stop("\n cannot use random effects with frequency >1")
	   }
	   # indicator for loss on capture 
	   write(model_data$imat$loc,con,ncolumns=n,append=TRUE)
	   write(t(model_data$time.intervals),con,ncolumns=nocc-1,append=TRUE)
	   write(ncol(model_data$Phi.dm),con,append=TRUE)
	   write(t(model_data$Phi.dm),con,ncolumns=ncol(model_data$Phi.dm),append=TRUE)
	   write(ncol(model_data$p.dm),con,append=TRUE)
	   write(t(model_data$p.dm),con,ncolumns=ncol(model_data$p.dm),append=TRUE)
	   if(model_data$Phi.fixed[1,1]== -1)
	   {
		   write(0,con,append=TRUE)
	   }else
	   {
		   index=(nocc-1)*(model_data$Phi.fixed[,1]-1)+model_data$Phi.fixed[,2]
		   write(nrow(model_data$Phi.fixed),con,append=TRUE)
		   write(rbind(index,model_data$Phi.fixed[,3]),con,ncolumns=2,append=TRUE)
	   }  
	   if(model_data$p.fixed[1,1]== -1)
	   {
		   write(0,con,append=TRUE)
	   }else
	   {
		   index=(nocc-1)*(model_data$p.fixed[,1]-1)+model_data$p.fixed[,2]
		   write(nrow(model_data$p.fixed),con,append=TRUE)
		   write(rbind(index,model_data$p.fixed[,3]),con,ncolumns=2,append=TRUE)
	   }  
	   close(con)
	   con=file(paste(tpl,".pin",sep=""),open="wt")
	   write(par$Phi,con,ncolumns=length(par$Phi),append=FALSE)
	   write(par$p,con,ncolumns=length(par$p),append=TRUE)
	   if(re) 
	   {
		   write(c(0.1,0.1),con,ncolumns=1,append=TRUE)
		   write(rep(0,n),con,ncolumns=n,append=TRUE)
		   write(rep(0,n),con,ncolumns=n,append=TRUE)
	   }
	   close(con)   
	   if(hessian)
	       xx=run_admb(tpl,extra.args=extra.args)
	   else
		   xx=run_admb(tpl,extra.args=paste(extra.args,"-nohess"))
	   convergence=attr(xx,"status")
	   if(is.null(convergence))convergence=0
	   res=read_admb(tpl)
	   cjs.beta.fixed=unscale.par(res$coefficients[1:(ncol(model_data$Phi.dm)+ncol(model_data$p.dm))],scale)
	   if(re)
	   {
		   cjs.beta.random=res$coefficients[(ncol(model_data$Phi.dm)+ncol(model_data$p.dm)+1):length(coef(res))]
		   names(cjs.beta.random)=paste("sigma_",names(cjs.beta.random),sep="")
	   }
	   else 
		   cjs.beta.random=NULL
	   cjs.beta=c(cjs.beta.fixed,cjs.beta.random)
	   beta=list(cjs.beta)
	   if(!is.null(res$hes))
	   {
		   beta.vcv=solvecov(res$hes)$inv
		   rownames(res$hes)=names(unlist(cjs.beta))
		   colnames(res$hes)=rownames(res$hes)
		   if(all(diag(beta.vcv>0))) 
		      res$cor=beta.vcv/outer(sqrt(diag(beta.vcv)),sqrt(diag(beta.vcv)))
	   }
	   else
		   beta.vcv=res$vcov
	   rownames(beta.vcv)=names(unlist(cjs.beta))
	   colnames(beta.vcv)=rownames(beta.vcv)
	   if(!re)
	   {
		   rownames(res$cor)=rownames(beta.vcv)
		   colnames(res$cor)=rownames(beta.vcv)
	   }
	   res$vcov=NULL
	   res=c(beta=beta,neg2lnl=-2*res$loglik,AIC=-2*res$loglik+2*res$npar,convergence=convergence,res)
	   res$beta.vcv=beta.vcv
   }
#  Restore non-accumulated, non-scaled dm's etc
   res$model_data=model_data.save
#  Assign S3 class values and return
   class(res)=c("crm","mle","cjs")
   return(res)
}


