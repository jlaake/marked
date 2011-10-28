cjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,Phi=NULL,p=NULL,initial=NULL,method,
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale, ...)
###################################################################################
# cjs - convenience function for optim to optimize the cjs likelihood (cjs.lnl) for
#       a given model specified by the Phi.dm and p.dm design matrices and
#       the data x.
#
# Arguments:
#    x               - processed dataframe created by process.data
#    ddl             - list of dataframes for design data
#    dml             - list of design matrices
#    parameters      - list of parameter model specifications
#    accumulate     - if TRUE will accumulate capture histories with common value
#                      and with a common design matrix and fixed values for Phi and p
#    Phi             - initial value for intercept 
#    p               - initial value for intercept
#    initial         - initial values for parameters if desired; if it is a
#                      named vector it will match with column names of design matrices
#    method          - method used in optim
#    hessian         - if TRUE, computes and returns hessian
#    debug           - if TRUE show iterations with par and -2lnl
#    chunk_size      - specifies amount of memory to use in accumulating capture histories
#                             use is 8*chunk_size/1e6 MB (default 80MB)
#    refit           - number of times to refit the model if it fails to converge
#    ...             - additional arguments passed to cjs.lnl or optim
#
# Value: list of results
#
###################################################################################
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
   model_data$Phi.dm=t(t(model_data$Phi.dm)/scale.phi)
   model_data$p.dm=t(t(model_data$p.dm)/scale.p)
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
	   mod=optimx(par,cjs.lnl,model_data=model_data,Phi.links=NULL,p.links=NULL,method=method,hessian=FALSE,
			         debug=debug,control=control,itnmax=itnmax,...)
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
	   res$vcv=cjs.hessian(res,Phi.links=NULL, p.links=NULL,all=FALSE)
   }   
#  Assign S3 class and return
   class(res)=c("crm","cjs")
   return(res)
}
cjs.hessian=function(model,Phi.links=NULL, p.links=NULL,all=FALSE)
{
	scale=c(model$scale$phi,model$scale$p)
	cat("\n Computing hessian\n")
	vcv=hessian(cjs.lnl,model$beta*scale,model_data=model$model_data,Phi.links=NULL, p.links=NULL,all=FALSE)
	vcv=try(solvecov(vcv))
	if(class(vcv)[1]=="try-error")
	{
		cat("\nUnable to invert hessian\n")
		return(NULL)
	}
	vcv=vcv$inv/outer(scale,scale,"*")
	colnames(vcv)=names(model$beta)
	rownames(vcv)=names(model$beta)
	return(vcv)
}

