#' Fitting function for Multistate CJS models with TMB
#' 
#' A function for computing MLEs for a Multi-state Cormack-Jolly-Seber open
#' population capture-recapture model for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{mscjs_tmb} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, S.dm,p.dm,Psi.dm,S.fixed,p.fixed,Psi.fixed and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for S and p to speed up execution
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
#' @param extra.args optional character string that is passed to admb 
#' @param clean if TRUE, deletes the tpl and executable files for amdb 
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
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Ford, J. H., M. V. Bravington, and J. Robbins. 2012. Incorporating individual variability into mark-recapture models. Methods in Ecology and Evolution 3:1047-1054.
mscjs_tmb=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		re=FALSE,compile=FALSE,extra.args="",clean=TRUE,...)
{
	accumulate=FALSE
	nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use x$intervals if none are in ddl$Phi
	if(!is.null(ddl$Phi$time.interval))
		time.intervals=matrix(ddl$Phi$time.interval,nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
	if(is.vector(x$time.intervals))
		time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
		time.intervals=x$time.intervals
#  Store data from x$data into x
	strata.labels=x$strata.labels
	uS=x$unobserved
	x=x$data
#  set default frequencies if not used
	freq=NULL
	if(!is.null(x$freq))freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
	ch=x$ch
	imat=process.ch(ch,freq,all=FALSE)
	chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=length(ch))
	for(nlabel in 1:length(strata.labels))
		chmat=t(apply(chmat,1,sub,pattern=strata.labels[nlabel],replacement=nlabel))
	chmat=t(apply(chmat,1,function(x) as.numeric(x)))
#  Use specified initial values or create if null
	if(is.null(initial))
		par=list(Psi=rep(0,ncol(dml$Psi$fe)),
				p=rep(0,ncol(dml$p$fe)),
				S=rep(0,ncol(dml$S$fe)))
	else
		par=set.initial(names(dml),dml,initial)$par
#  Create list of model data for optimization
	model_data=list(S.dm=dml$S$fe,p.dm=dml$p$fe,Psi.dm=dml$Psi$fe,imat=imat,S.fixed=parameters$S$fixed,
			p.fixed=parameters$p$fixed,Psi.fixed=parameters$Psi$fixed,time.intervals=time.intervals)
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
	scale=set.scale(names(dml),model_data,scale)
	model_data=scale.dm(model_data,scale)
	setup_tmb("multistate_tmb",clean=clean)
	cat("\nbuilding TMB program\n")

	# S design matrix
	phidm=as.matrix(model_data$S.dm)
	phifix=rep(-1,nrow(phidm))
	if(!is.null(ddl$S$fix))
		phifix[!is.na(ddl$S$fix)]=ddl$S$fix[!is.na(ddl$S$fix)]
	phi_slist=simplify_indices(cbind(phidm,phifix))
	
	# p design matrix
	pdm=as.matrix(model_data$p.dm)
	pfix=rep(-1,nrow(pdm))
	if(!is.null(ddl$p$fix))
		pfix[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
	p_slist=simplify_indices(cbind(pdm,pfix))

	#Psi design matrix
	psidm=as.matrix(model_data$Psi.dm)
	psifix=rep(-1,nrow(psidm))
	if(!is.null(ddl$Psi$fix))
		psifix[!is.na(ddl$Psi$fix)]=ddl$Psi$fix[!is.na(ddl$Psi$fix)]
	psi_slist=simplify_indices(cbind(psidm,psifix))
	
	f = MakeADFun(data=list(n=length(model_data$imat$freq),m=model_data$imat$nocc,nS=length(strata.labels),
					ch=chmat,frst=model_data$imat$first,freq=model_data$imat$freq,tint=model_data$time.intervals,
					kphi=ncol(phidm),nrowphi=length(phi_slist$set),	phidm=phidm[phi_slist$set,,drop=FALSE],
					phifix=phifix[phi_slist$set],phiindex=phi_slist$indices[ddl$S.indices],
					kp=ncol(pdm),nrowp=length(p_slist$set),pdm=pdm[p_slist$set,,drop=FALSE],
					pfix=pfix[p_slist$set],pindex=p_slist$indices[ddl$p.indices],
					kpsi=ncol(psidm),nrowpsi=length(psi_slist$set),	psidm=psidm[psi_slist$set,,drop=FALSE],
					psifix=psifix[psi_slist$set],psiindex=psi_slist$indices[ddl$Psi.indices]),
			        parameters=list(phibeta=par$S,pbeta=par$p,psibeta=par$Psi)
					,DLL="multistate_tmb")
	cat("\nrunning TMB program\n")                         
	if(method=="nlminb")
	{
		mod=nlminb(f$par,f$fn,f$gr,control=control,itnmax=itnmax,...)
		lnl=mod$objective
		par=mod$par
		convergence=mod$convergence
	} else
	{
		control$starttests=FALSE
		mod=optimx(f$par,f$fn,f$gr,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
		par <- coef(mod, order="value")[1, ]
		mod=as.list(summary(mod, order="value")[1, ])
	    convergence=mod$convcode
		lnl=mod$value		
	}
	fixed.npar=(ncol(phidm)+ncol(pdm)-2)
			if(p_nre+phi_nre>0)
			{
				if(getreals) 
					par_summary=sdreport(f,getReportCovariance=FALSE)
				else
					par_summary=sdreport(f,getJointPrecision=TRUE)
				
				par=par_summary$par.fixed[1:fixed.npar]
				cjs.beta.fixed=unscale.par(par,scale)
				cjs.beta.sigma=par_summary$par.fixed[-(1:fixed.npar)]
				sigma=NULL
				if(phi_krand>0)
				{
					Phi_sigma=cjs.beta.sigma[1:phi_krand]
					names(Phi_sigma)=colnames(phi_randDM)
					sigma=list(Phi_logsigma=Phi_sigma)
				} 
				if(p_krand>0)
				{
					p_sigma=cjs.beta.sigma[(phi_krand+1):(phi_krand+p_krand)]
					names(p_sigma)=colnames(p_randDM)
					sigma=c(sigma,list(p_logsigma=p_sigma))
				}
				cjs.beta=c(cjs.beta.fixed,sigma)
				beta.vcv=par_summary$cov.fixed
			}else
			{	
				cjs.beta=unscale.par(par,scale)
				if(hessian) 
				{
					message("Computing hessian...")
					beta.vcv=solve(f$he(par))
					colnames(beta.vcv)=names(unlist(cjs.beta))
					rownames(beta.vcv)=colnames(beta.vcv)
				} else
					beta.vcv=NULL
			}	
			# Create results list 
			if(getreals&p_nre+phi_nre>0)
			{
				reals=split(par_summary$value,names(par_summary$value))
				reals.se=split(par_summary$sd,names(par_summary$value))	
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
			
#  Restore non-accumulated, non-scaled dm's etc
			res$model_data=model_data.save
#  Assign S3 class values and return
			class(res)=c("crm","admb","mle","mscjs")
			return(res)
			}




