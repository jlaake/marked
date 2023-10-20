#' Fitting function for Multistate CJS live-dead models with TMB
#' 
#' A function for computing MLEs for a Multi-state Cormack-Jolly-Seber open
#' population capture-recapture with dead recoveries for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{msld_tmb} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of simplified dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, S.dm,r.dm,p.dm,Psi.dm,S.fixed,r.fixed,p.fixed,Psi.fixed and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for all parameters speed up execution
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
#' @param re if TRUE creates random effect model admbcjsre.tpl and runs admb optimizer
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to tmb 
#' @param clean if TRUE, deletes the dll and recompiles 
#' @param getreals if TRUE, compute real values and std errors for TMB models; may want to set as FALSE until model selection is complete
#' @param useHess if TRUE, the TMB hessian function is used for optimization; using hessian is typically slower with many parameters but can result in a better solution
#' @param savef if TRUE, save optimization function in model for reporting
#' @param ... not currently used
#' @export
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optim}; if 0 \code{optim} thinks it
#' converged} \item{count}{\code{optim} results of number of function
#' evaluations} \item{reals}{dataframe of data and real S and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laak
smsld_tmb=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		re=FALSE,compile=FALSE,extra.args="",clean=FALSE,getreals=FALSE, useHess=FALSE,savef=FALSE,...)
{
  # load fullddl
  fullddl=NULL
  load("tmp.rda")
	accumulate=FALSE
	nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use x$intervals if none are in ddl$Phi
	if(!is.null(ddl$S$time.interval))		
	  time.intervals=matrix(fullddl$S$time.interval[fullddl$S$stratum==x$strata.labels[1]],nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
	if(is.vector(x$time.intervals))
		time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
		time.intervals=x$time.intervals

#  Store data from x$data into x
	strata.labels=c(x$strata.labels,"1")
	uS=x$unobserved
	x=x$data
#  set default frequencies if not used
	freq=NULL
	if(!is.null(x$freq))freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
	ch=x$ch
	imat=process.ch(ch,freq,all=FALSE)
	chmat=matrix((unlist(strsplit(ch,""))),byrow=TRUE,ncol=nocc,nrow=length(ch))
	for(i in 1:length(strata.labels))
	{
	  nlabel=length(strata.labels)-i+1
	  chmat=t(apply(chmat,1,sub,pattern=strata.labels[nlabel],replacement=nlabel))
	}
	chmat=t(apply(chmat,1,function(x) as.numeric(x)))
#  Use specified initial values or create if null
	if(is.null(initial))
		par=list(Psi=rep(0,ncol(dml$Psi$fe)),
				p=rep(0,ncol(dml$p$fe)),
				r=rep(0,ncol(dml$r$fe)),
				S=rep(0,ncol(dml$S$fe)))
	else
		par=set.initial(names(dml),dml,initial)$par
#  Create list of model data for optimization
	model_data=list(S.dm=dml$S$fe,r.dm=dml$r$fe,p.dm=dml$p$fe,Psi.dm=dml$Psi$fe,imat=imat,S.fixed=parameters$S$fixed,
	                r.fixed=parameters$r$fixed,p.fixed=parameters$p$fixed,Psi.fixed=parameters$Psi$fixed,
	                time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
	if(accumulate)
	{
		cat("Accumulating capture histories based on design. This can take awhile.\n")
		flush.console()
		model_data.save=model_data   
		#model_data=mscjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   S.links=create.links(S.dm)
#   S.links=which(S.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#  Scale the design matrices and parameters with either input scale or computed scale
	scale=1
	scale=set_scale(names(dml),model_data,scale)
	model_data=scale_dm(model_data,scale)
	setup_tmb("smsld_tmb",clean=clean)

	# S design matrix
	phidm=as.matrix(model_data$S.dm)
	phifix=rep(-1,nrow(phidm))
	if(!is.null(ddl$S$fix))
		phifix[!is.na(ddl$S$fix)]=ddl$S$fix[!is.na(ddl$S$fix)]
	phi_slist=simplify_indices(cbind(phidm,phifix))
	phimixed=mixed.model.admb(parameters$S$formula,fullddl$S)
	nphisigma=0
	if(!is.null(phimixed$re.dm))nphisigma=ncol(phimixed$re.dm)
	if(!is.null(phimixed$re.dm))
	{
		phimixed$re.indices[fullddl$S$Time<fullddl$S$Cohort,]=NA
		phimixed=reindex(phimixed,fullddl$S$id)
		
		# random effect data
		phi_krand=ncol(phimixed$re.dm)
		phi_randDM=phimixed$re.dm
		phi_randIndex=phimixed$re.indices
		phi_counts=phimixed$index.counts
		mx=max(phimixed$index.counts)
		phi_idIndex=t(sapply(phimixed$used.indices,function(x) return(c(x,rep(0,mx-length(x))))))
		phi_nre=max(phi_idIndex)
		if(nrow(phi_idIndex)==1)phi_idIndex=t(phi_idIndex)
		if(phi_krand==1)phi_randIndex=matrix(as.vector(t(phi_randIndex)),ncol=1)
#   simplify random indices and dm
		s1=simplify_indices(phi_idIndex)
		phi_idIndex=phi_idIndex[s1$set,]
		phi_idIndex_i=s1$indices
		s1=simplify_indices(phi_randIndex)
		phi_randIndex=phi_randIndex[s1$set,]
		phi_randIndex_i=s1$indices
		s1=simplify_indices(phi_randDM)
		phi_randDM=phi_randDM[s1$set,]
		phi_randDM_i=s1$indices
	} else {
		phi_nre=0
		phi_krand=0
		phi_randDM=matrix(0,nrow=0,ncol=0)
		phi_randIndex=matrix(0,nrow=0,ncol=0)
		phi_counts=vector("integer",length=0)
		phi_idIndex=matrix(0,nrow=0,ncol=0)
	}

	# r design matrix
	rdm=as.matrix(model_data$r.dm)
	rfix=rep(-1,nrow(rdm))
	if(!is.null(ddl$r$fix))
	  rfix[!is.na(ddl$r$fix)]=ddl$r$fix[!is.na(ddl$r$fix)]
	r_slist=simplify_indices(cbind(rdm,rfix))
	rmixed=mixed.model.admb(parameters$r$formula,fullddl$r)
	nrsigma=0
	if(!is.null(rmixed$re.dm))nrsigma=ncol(rmixed$re.dm)
	if(!is.null(rmixed$re.dm))
	{
	  rmixed$re.indices[fullddl$r$Time<fullddl$r$Cohort,]=NA
	  rmixed=reindex(rmixed,fullddl$r$id)
	  
	  # random effect data
	  r_krand=ncol(rmixed$re.dm)
	  r_randDM=rmixed$re.dm
	  r_randIndex=rmixed$re.indices
	  r_counts=rmixed$index.counts
	  mx=max(rmixed$index.counts)
	  r_idIndex=t(sapply(rmixed$used.indices,function(x) return(c(x,rep(0,mx-length(x))))))
	  r_nre=max(r_idIndex)
	  if(nrow(r_idIndex)==1)r_idIndex=t(r_idIndex)
	  if(r_krand==1)r_randIndex=matrix(as.vector(t(r_randIndex)),ncol=1)
	} else {
	  r_nre=0
	  r_krand=0
	  r_randDM=matrix(0,nrow=0,ncol=0)
	  r_randIndex=matrix(0,nrow=0,ncol=0)
	  r_counts=vector("integer",length=0)
	  r_idIndex=matrix(0,nrow=0,ncol=0)
	}
	
	
	# p design matrix
	pdm=as.matrix(model_data$p.dm)
	pfix=rep(-1,nrow(pdm))
	if(!is.null(ddl$p$fix))
		pfix[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
	p_slist=simplify_indices(cbind(pdm,pfix))
	pmixed=mixed.model.admb(parameters$p$formula,fullddl$p)
	npsigma=0
	if(!is.null(pmixed$re.dm))npsigma=ncol(pmixed$re.dm)
	
	if(!is.null(pmixed$re.dm))
	{
		pmixed$re.indices[fullddl$p$Time<fullddl$Cohort,]=NA
		pmixed=reindex(pmixed,fullddl$p$id)
		# random effect data	
		p_krand=ncol(pmixed$re.dm)
		p_randDM=pmixed$re.dm
		p_randIndex=pmixed$re.indices
		p_counts=pmixed$index.counts
		mx=max(pmixed$index.counts)
		p_idIndex=t(sapply(pmixed$used.indices,function(x) return(c(x,rep(0,mx-length(x))))))
		if(nrow(p_idIndex)==1)p_idIndex=t(p_idIndex)
		p_nre=max(p_idIndex)
		if(p_krand==1)p_randIndex=matrix(as.vector(t(p_randIndex)),ncol=1)
	} else {
		p_nre=0
		p_krand=0
		p_randDM=matrix(0,nrow=0,ncol=0)
		p_randIndex=matrix(0,nrow=0,ncol=0)
		p_counts=vector("integer",length=0)
		p_idIndex=matrix(0,nrow=0,ncol=0)
	}
	
	#Psi design matrix
	psidm=as.matrix(model_data$Psi.dm)
	psifix=rep(-1,nrow(psidm))
	if(!is.null(ddl$Psi$fix))
		psifix[!is.na(ddl$Psi$fix)]=ddl$Psi$fix[!is.na(ddl$Psi$fix)]
	psi_slist=simplify_indices(cbind(psidm,psifix))
	psimixed=mixed.model.admb(parameters$Psi$formula,fullddl$Psi)
	npsisigma=0
	if(!is.null(psimixed$re.dm))npsisigma=ncol(psimixed$re.dm)
	
	if(!is.null(psimixed$re.dm))
	{
		psimixed$re.indices[fullddl$Psi$Time<fullddl$Cohort,]=NA
		psimixed=reindex(psimixed,fullddl$Psi$id)
		# random effect data	
		psi_krand=ncol(psimixed$re.dm)
		psi_randDM=psimixed$re.dm
		psi_randIndex=psimixed$re.indices
		psi_counts=psimixed$index.counts
		mx=max(psimixed$index.counts)
		psi_idIndex=t(sapply(psimixed$used.indices,function(x) return(c(x,rep(0,mx-length(x))))))
		if(nrow(psi_idIndex)==1)psi_idIndex=t(psi_idIndex)
		psi_nre=max(psi_idIndex)
		if(psi_krand==1)psi_randIndex=matrix(as.vector(t(psi_randIndex)),ncol=1)
	} else {
		psi_nre=0
		psi_krand=0
		psi_randDM=matrix(0,nrow=0,ncol=0)
		psi_randIndex=matrix(0,nrow=0,ncol=0)
		psi_counts=vector("integer",length=0)
		psi_idIndex=matrix(0,nrow=0,ncol=0)
	}
	rm(fullddl)
	gc()
	message("Building optimization function")
	f = MakeADFun(data=list(n=length(model_data$imat$freq),m=model_data$imat$nocc,nS=length(strata.labels)-1,
					ch=chmat,frst=model_data$imat$first,freq=model_data$imat$freq,tint=model_data$time.intervals,
					nrowphi=length(phi_slist$set),	phidm=phidm[phi_slist$set,,drop=FALSE],
					phifix=phifix[phi_slist$set],phiindex=phi_slist$indices[ddl$S.indices],
					phi_nre=phi_nre,phi_krand=phi_krand,phi_randDM=phi_randDM,phi_randDM_i=phi_randDM_i,
					phi_randIndex=phi_randIndex,phi_randIndex_i=phi_randIndex_i,phi_counts=phi_counts,phi_idIndex=phi_idIndex,
					phi_idIndex_i=phi_idIndex_i,nrowr=length(r_slist$set),	rdm=rdm[r_slist$set,,drop=FALSE],
					rfix=rfix[r_slist$set],rindex=r_slist$indices[ddl$r.indices],
					r_nre=r_nre,r_krand=r_krand,r_randDM=r_randDM,
					r_randIndex=r_randIndex,r_counts=r_counts,r_idIndex=r_idIndex,
					nrowp=length(p_slist$set),pdm=pdm[p_slist$set,,drop=FALSE],
					pfix=pfix[p_slist$set],pindex=p_slist$indices[ddl$p.indices],
					p_nre=p_nre,p_krand=p_krand,p_randDM=p_randDM,
					p_randIndex=p_randIndex,p_counts=p_counts,p_idIndex=p_idIndex,
					nrowpsi=length(psi_slist$set),	psidm=psidm[psi_slist$set,,drop=FALSE],
					psifix=psifix[psi_slist$set],psiindex=psi_slist$indices[ddl$Psi.indices],
					psi_nre=psi_nre,psi_krand=psi_krand,psi_randDM=psi_randDM,
					psi_randIndex=psi_randIndex,psi_counts=psi_counts,psi_idIndex=psi_idIndex,getreals=as.integer(getreals)),
			        parameters=list(phibeta=par$S,rbeta=par$r,pbeta=par$p,psibeta=par$Psi,log_sigma_phi=rep(-1,nphisigma),
			                   log_sigma_r=rep(-1,nrsigma),log_sigma_p=rep(-1,npsigma),log_sigma_psi=rep(-1,npsisigma),
			                   u_phi=rep(0,phi_nre),u_r=rep(0,r_nre), u_p=rep(0,p_nre),u_psi=rep(0,psi_nre)),
					                random=c("u_phi","u_r","u_p","u_psi"),DLL="smsld_tmb")
	cat("\nrunning TMB program\n")                         
	if(method=="nlminb")
	{
	  if(!useHess)
		   mod=nlminb(f$par,f$fn,f$gr,control=control,...)
	  else
	    mod=nlminb(f$par,f$fn,f$gr,f$he,control=control,...)
		lnl=mod$objective
		par=mod$par
		convergence=mod$convergence
	} else
	{
		if(method=="SANN")
		{		  
		  control$maxit=itnmax
		  mod=optim(f$par,f$fn,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
		  par=mod$par
		  convergence=mod$convergence
    } else
    {
      control$starttests=FALSE
      if(!useHess)
        mod=optimx(f$par,f$fn,f$gr,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
      else
        mod=optimx(f$par,f$fn,f$gr,f$he,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
      par <- coef(mod, order="value")[1, ]
      mod=as.list(summary(mod, order="value")[1, ])
      convergence=mod$convcode
    }
		lnl=mod$value		
	}
	fixed.npar=ncol(phidm)+ncol(rdm)+ncol(pdm)+ncol(psidm)
	load("tmp.rda")
	if(getreals)
	  par_summary=sdreport(f,getReportCovariance=FALSE)
	else
	  par_summary=sdreport(f,getJointPrecision=TRUE)
	if(p_nre+phi_nre>0)
	{
	 	par=par_summary$par.fixed[1:fixed.npar]
		cjs.beta.fixed=unscale_par(par,scale)
		cjs.beta.sigma=par_summary$par.fixed[-(1:fixed.npar)]
		sigma=NULL
		if(phi_krand>0)
		{
			Phi_sigma=cjs.beta.sigma[1:phi_krand]
			names(Phi_sigma)=colnames(phi_randDM)
			sigma=list(Phi_logsigma=Phi_sigma)
		} 
		if(r_krand>0)
		{
		  r_sigma=cjs.beta.sigma[1:r_krand]
		  names(r_sigma)=colnames(r_randDM)
		  sigma=list(r_logsigma=r_sigma)
		} 
		if(p_krand>0)
		{
			p_sigma=cjs.beta.sigma[(phi_krand+1):(phi_krand+p_krand)]
			names(p_sigma)=colnames(p_randDM)
			sigma=c(sigma,list(p_logsigma=p_sigma))
		}
		if(psi_krand>0)
		{
			Psi_sigma=cjs.beta.sigma[(phi_krand+p_krand+1):(phi_krand+p_krand+psi_krand)]
			names(Psi_sigma)=colnames(psi_randDM)
			sigma=c(sigma,list(Psi_logsigma=Psi_sigma))
		} 
		cjs.beta=c(cjs.beta.fixed,sigma)
		beta.vcv=par_summary$cov.fixed
	}else
	{	
		cjs.beta=unscale_par(par,scale)
		if(hessian) 
		{
			message("Computing hessian...")
			beta.vcv=solvecov(f$he(par))$inv
			colnames(beta.vcv)=names(unlist(cjs.beta))
			rownames(beta.vcv)=colnames(beta.vcv)
		} else
			beta.vcv=NULL
	}	
	# Create results list 
	if(getreals)
	{
		reals=split(par_summary$value,names(par_summary$value))
		reals.se=split(par_summary$sd,names(par_summary$value))	
		names(reals)=c("p","S","Psi")
		names(reals.se)=c("p","S","Psi")
		reals$S[fullddl$S$Time<fullddl$S$Cohort]=NA
		reals.se$S[fullddl$S$Time<fullddl$S$Cohort]=NA
		reals$p[fullddl$p$Time<fullddl$p$Cohort]=NA
		reals.se$p[fullddl$p$Time<fullddl$p$Cohort]=NA
		reals$Psi=as.vector(aperm(array(reals$Psi,dim=c(model_data$imat$nocc-1,length(strata.labels),length(strata.labels),length(model_data$imat$freq))),c(3,2,1,4)))
		reals.se$Psi=as.vector(aperm(array(reals.se$Psi,dim=c(model_data$imat$nocc-1,length(strata.labels),length(strata.labels),length(model_data$imat$freq))),c(3,2,1,4)))
#		reals$Psi=as.vector(aperm(array(reals$Psi,dim=c(model_data$imat$nocc-1,(length(strata.labels))^2,length(model_data$imat$freq))),c(2,1,3)))
#		reals.se$Psi=as.vector(aperm(array(reals.se$Psi,dim=c(model_data$imat$nocc-1,(length(strata.labels))^2,length(model_data$imat$freq))),c(2,1,3)))
		reals$Psi[fullddl$Psi$Time<fullddl$Psi$Cohort]=NA
		reals.se$Psi[fullddl$Psi$Time<fullddl$Psi$Cohort]=NA
	}
	else
	{
		reals=NULL
		reals.se=NULL
	}
	res=list(beta=cjs.beta,neg2lnl=2*lnl,AIC=2*lnl+2*sum(sapply(cjs.beta,length)),
			beta.vcv=beta.vcv,reals=reals,reals.se=reals.se,convergence=convergence,optim.details=mod,
			model_data=model_data,
			options=list(scale=scale,accumulate=accumulate,initial=initial,method=method,
			chunk_size=chunk_size,itnmax=itnmax,control=control))		
		
# Restore non-accumulated, non-scaled dm's etc
	res$model_data=model_data.save
# if savef add it to the model
	if(savef)res$f=f
# Assign S3 class values and return
	class(res)=c("crm","admb","mle","mscjs")
	return(res)
}




