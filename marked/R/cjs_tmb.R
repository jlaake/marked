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
#' @param crossed if TRUE it uses cjs.tpl or cjs_reml.tpl if reml=FALSE or TRUE respectively; if FALSE, then it uses cjsre which can use Gauss-Hermite integration
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param reml if set to TRUE uses cjs_reml if crossed 
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param getreals  if TRUE, will compute real Phi and p values and std errors
#' @param prior if TRUE will expect vectors of prior values in list prior.list
#' @param prior.list which contains for normal distributions 1) mu_phi_prior: vector of mu values for phi_beta, 2) sigma_phi_prior: vector of sigma values for phi_beta,
#' 3) mu_p_prior: vector of mu values for p_beta, 4) sigma_p_prior: vector of sigma values for p_beta, 5) random_mu_phi_prior: vector of mu values for ln sigma of random effects, 
#' 6) random_sigma_phi_prior: vector of sigma values for ln sigma_phi, 7) random_mu_p_prior: vector of mu values for ln sigma_p, 8) random_sigma_p_prior: vector of sigma values for ln sigma_p. 
#' @param tmbfct either "f1" - default or "f2" - any random effects treated as fixed effects or "f3" fixed effects fixed at mode and no random effects.
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optim} or \code{\link{cjs.lnl}}
#' @import R2admb optimx TMB
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' Elements are 1) beta: named vector of parameter estimatesm 2) lnl: -2*log
#' likelihood, 3) AIC: lnl + 2* number of parameters, 4) convergence: result from \code{optim}; if 0 \code{optim} thinks it
#' converged, 5) count:\code{optim} results of number of function
#' evaluations, 6) reals: dataframe of data and real Phi and p estimates for
#' each animal-occasion excluding those that occurred before release, 7) vcv:var-cov matrix of betas if hessian=TRUE was set.
#' @author Jeff Laake 
#' @references Pledger, S., K. H. Pollock, et al. (2003). Open
#' capture-recapture models with heterogeneity: I. Cormack-Jolly-Seber model.
#' Biometrics 59(4):786-794.
cjs_tmb=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		crossed=TRUE,compile=TRUE,extra.args=NULL,reml,clean=FALSE,getreals=FALSE,prior=FALSE,
		prior.list=NULL,tmbfct="f1",...)
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
#  Create fixed matrices in parameters
	parameters=create.fixed.matrix(ddl,parameters)
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
		par=set.initial(names(dml),dml,initial)$par
	initial=par
#  Create list of model data for optimization
	model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe,imat=imat,Phi.fixed=parameters$Phi$fixed,
			p.fixed=parameters$p$fixed,time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
#   Problems with accumulation and fixed values 10 Jan 2014; turned off accumulate if fixed
	if(parameters$p$fixed[1,1]>0 | parameters$Phi$fixed[1,1]>0) accumulate=FALSE
	if(accumulate)
	{
		message("Accumulating capture histories based on design. This can take awhile...")
		model_data.save=model_data   
		model_data=cjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#   Scale the design matrices and parameters with either input scale or computed scale of mean
#   Currently no scaling
    scale=1
	scale=set.scale(names(dml),model_data,scale)
	model_data=scale.dm(model_data,scale)

		########################################################################
#      CJS with TMB
		########################################################################
		#phi dm portion
 		phidm=model_data$Phi.dm
		phidm=cbind(phidm,rep(-1,nrow(phidm)))
		if(model_data$Phi.fixed[1,1]!= -1)
		{
			index=(nocc-1)*(model_data$Phi.fixed[,1]-1)+model_data$Phi.fixed[,2]
			phidm[index,ncol(phidm)]=model_data$Phi.fixed[,3]
			phidm[index,1:(ncol(phidm)-1)]=0
		}
		phimixed=mixed.model.admb(parameters$Phi$formula,ddl$Phi)
		nphisigma=0
		if(!is.null(phimixed$re.dm))nphisigma=ncol(phimixed$re.dm)
		if(!is.null(phimixed$re.dm))
		{
			phimixed$re.indices[ddl$Phi$Time<ddl$Phi$Cohort,]=NA
			phimixed=reindex(phimixed,ddl$Phi$id)

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
			phi_refreq=phimixed$freq
		} else {
			phi_nre=0
			phi_krand=0
			phi_randDM=matrix(0,nrow=0,ncol=0)
			phi_randIndex=matrix(0,nrow=0,ncol=0)
			phi_counts=vector("integer",length=0)
			phi_idIndex=matrix(0,nrow=0,ncol=0)
			phi_refreq=vector("numeric",length=0)
		}
		
		# p dm portion        
        mlist=proc.form(parameters$p$formula)
		p_re_names=NULL
        if(!is.null(mlist$re.model))
	       p_re_names=sub("^\\s+", "",sapply(strsplit(names(mlist$re.model),"\\|"),function(x)x[2]))
		pdm=model_data$p.dm
		pdm=cbind(pdm,rep(-1,nrow(pdm)))
		if(model_data$p.fixed[1,1]!= -1)
		{
			index=(nocc-1)*(model_data$p.fixed[,1]-1)+model_data$p.fixed[,2]-1
			pdm[index,ncol(pdm)]=model_data$p.fixed[,3]
			pdm[index,1:(ncol(pdm)-1)]=0
		}
		pmixed=mixed.model.admb(parameters$p$formula,ddl$p)
		npsigma=0
		if(!is.null(pmixed$re.dm))npsigma=ncol(pmixed$re.dm)
		
		if(!is.null(pmixed$re.dm))
		{
			pmixed$re.indices[ddl$p$Time<ddl$Cohort,]=NA
			pmixed=reindex(pmixed,ddl$p$id)
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
			p_refreq=pmixed$freq
		} else {
			p_nre=0
			p_krand=0
			p_randDM=matrix(0,nrow=0,ncol=0)
			p_randIndex=matrix(0,nrow=0,ncol=0)
			p_counts=vector("integer",length=0)
			p_idIndex=matrix(0,nrow=0,ncol=0)
			p_refreq=vector("numeric",length=0)
		}
	
		# if priors are desired, set up appropriate structure
		if(prior)
		{
        # priors for fixed effects
			if(is.null(prior.list)|| is.null(prior.list$mu_phi_prior)) 
				mu_phi_prior=0
			else
				mu_phi_prior=prior.list$mu_phi_prior
		  	if(length(mu_phi_prior)==1) mu_phi_prior=rep(mu_phi_prior,ncol(phidm))
		   	if(length(mu_phi_prior)!=ncol(phidm))stop("\nMismatch between length of mu_phi_prior and number of phi_betas")
		   	if(is.null(prior.list)|| is.null(prior.list$mu_p_prior)) 
				mu_p_prior=0
			else
				mu_p_prior=prior.list$mu_p_prior
		   	if(length(mu_p_prior)==1) mu_p_prior=rep(mu_p_prior,ncol(pdm))
		   	if(length(mu_p_prior)!=ncol(pdm))stop("\nMismatch between length of mu_p_prior and number of p_betas")
		   	if(is.null(prior.list)|| is.null(prior.list$sigma_phi_prior)) 
				sigma_phi_prior=100
			else
				sigma_phi_prior=prior.list$sigma_phi_prior
		   	if(length(sigma_phi_prior)==1) sigma_phi_prior=rep(sigma_phi_prior,ncol(phidm))
		   	if(length(sigma_phi_prior)!=ncol(phidm))stop("\nMismatch between length of sigma_phi_prior and number of phi_betas")
		   	if(is.null(prior.list)|| is.null(prior.list$sigma_p_prior)) 
				sigma_p_prior=100
			else
				sigma_p_prior=prior.list$sigma_p_prior
		   	if(length(sigma_p_prior)==1) sigma_p_prior=rep(sigma_p_prior,ncol(pdm))
		   	if(length(sigma_p_prior)!=ncol(pdm))stop("\nMismatch between length of sigma_p_prior and number of p_betas")
		# priors for random effects  
	        if(phi_krand>0)
			{	
				if(is.null(prior.list)|| is.null(prior.list$random_mu_phi_prior)) 
					random_mu_phi_prior=0
				else
					random_mu_phi_prior=prior.list$random_mu_phi_prior
				if(length(random_mu_phi_prior)==1) random_mu_phi_prior=rep(random_mu_phi_prior,phi_krand)
				if(length(random_mu_phi_prior)!=phi_krand)stop("\nMismatch between length of random_mu_phi_prior and number of phi random effects")
				if(is.null(prior.list)|| is.null(prior.list$random_sigma_phi_prior)) 
					random_sigma_phi_prior=1
				else
					random_sigma_phi_prior=prior.list$random_sigma_phi_prior
				if(length(random_sigma_phi_prior)==1) random_sigma_phi_prior=rep(random_sigma_phi_prior,phi_krand)
				if(length(random_sigma_phi_prior)!=phi_krand)stop("\nMismatch between length of random_sigma_phi_prior and number of phi random effects")
			} else {
				random_mu_phi_prior=vector("numeric",length=0)
				random_sigma_phi_prior=vector("numeric",length=0)
			}
			if(phi_krand>0)
			{	
				if(is.null(prior.list)|| is.null(prior.list$random_mu_p_prior)) 
					random_mu_p_prior=0
				else
					random_mu_p_prior=prior.list$random_mu_p_prior
			 	if(length(random_mu_p_prior)==1) random_mu_p_prior=rep(random_mu_p_prior,p_krand)
				if(length(random_mu_p_prior)!=p_krand)stop("\nMismatch between length of random_mu_p_prior and number of p random effects")
				if(is.null(prior.list)|| is.null(prior.list$random_sigma_p_prior)) 
					random_sigma_p_prior=1
				else
					random_sigma_p_prior=prior.list$random_sigma_p_prior
				if(length(random_sigma_p_prior)==1) random_sigma_p_prior=rep(random_sigma_p_prior,p_krand)
				if(length(random_sigma_p_prior)!=p_krand)stop("\nMismatch between length of random_sigma_p_prior and number of p random effects")
			} else{
				random_mu_p_prior=vector("numeric",length=0)
				random_sigma_p_prior=vector("numeric",length=0)
			}
		} else {
		# no priors but need to create empty vectors
	       mu_phi_prior=vector("numeric",length=0)
		   sigma_phi_prior=vector("numeric",length=0)
		   mu_p_prior=vector("numeric",length=0)
		   sigma_p_prior=vector("numeric",length=0)
		   random_mu_phi_prior=vector("numeric",length=0)
		   random_sigma_phi_prior=vector("numeric",length=0)
		   random_mu_p_prior=vector("numeric",length=0)
		   random_sigma_p_prior=vector("numeric",length=0)
	   }	
		setup_tmb("cjsre_tmb",clean=clean)
		cat("\nbuilding TMB program\n")                         
		# Create AD function with data and parameters
        # With INLA type approach will need to run MakeADFun 3x. 
        # f1 - function with random= u's - optimize
        # f2 - function with u's not random - optimize
        # f3 - function with fixed effect parameter values specified at MLEs; use MAP to fix parameters; us not random - optimize
        # https://github.com/James-Thorson/2016_Spatio-temporal_models/issues/8
		if(tmbfct=="f1")
		{
			f = MakeADFun(data=list(m=model_data$imat$nocc,ch=model_data$imat$chmat,freq=model_data$imat$freq,frst=model_data$imat$first,
							lst=model_data$imat$last,loc=model_data$imat$loc,tint=model_data$time.intervals,
							phi_fixedDM=phidm,phi_nre=phi_nre,phi_krand=phi_krand,phi_randDM=phi_randDM,
							phi_randIndex=phi_randIndex,phi_counts=phi_counts,phi_idIndex=phi_idIndex,
							p_fixedDM=pdm,p_nre=p_nre,p_krand=p_krand,p_randDM=p_randDM,
							p_randIndex=p_randIndex,p_counts=p_counts,p_idIndex=p_idIndex,getreals=as.integer(getreals),
							prior=as.numeric(prior),mu_phi_prior=mu_phi_prior,sigma_phi_prior=sigma_phi_prior,
							mu_p_prior=mu_p_prior,sigma_p_prior=sigma_p_prior,random_mu_phi_prior=random_mu_phi_prior,
							random_sigma_phi_prior=random_sigma_phi_prior,random_mu_p_prior=mu_p_prior,random_sigma_p_prior=sigma_p_prior),
					        parameters=list(phi_beta=initial$Phi,p_beta=initial$p,
							log_sigma_phi=rep(-1,nphisigma),log_sigma_p=rep(-1,npsigma),u_phi=rep(0,phi_nre),u_p=rep(0,p_nre)),
					        random=c("u_phi","u_p"),DLL="cjsre_tmb")
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
		
#               Restore non-accumulated, non-scaled dm's etc
	            res$model_data=model_data.save
				res$tmbfct=f
#               Assign S3 class values and return
	            class(res)=c("crm","mle","cjs")
	            return(res)
		}
		if(tmbfct=="f2")
		{
			f = MakeADFun(data=list(m=model_data$imat$nocc,ch=model_data$imat$chmat,freq=model_data$imat$freq,frst=model_data$imat$first,
							lst=model_data$imat$last,loc=model_data$imat$loc,tint=model_data$time.intervals,
							phi_fixedDM=phidm,phi_nre=phi_nre,phi_krand=phi_krand,phi_randDM=phi_randDM,
							phi_randIndex=phi_randIndex,phi_counts=phi_counts,phi_idIndex=phi_idIndex,
							p_fixedDM=pdm,p_nre=p_nre,p_krand=p_krand,p_randDM=p_randDM,
							p_randIndex=p_randIndex,p_counts=p_counts,p_idIndex=p_idIndex,getreals=as.integer(getreals),
							prior=as.numeric(prior),mu_phi_prior=mu_phi_prior,sigma_phi_prior=sigma_phi_prior,
							mu_p_prior=mu_p_prior,sigma_p_prior=sigma_p_prior,random_mu_phi_prior=random_mu_phi_prior,
							random_sigma_phi_prior=random_sigma_phi_prior,random_mu_p_prior=mu_p_prior,random_sigma_p_prior=sigma_p_prior),
					parameters=list(phi_beta=initial$Phi,p_beta=initial$p,
							log_sigma_phi=rep(-1,nphisigma),log_sigma_p=rep(-1,npsigma),u_phi=rep(0,phi_nre),u_p=rep(0,p_nre)),
					,DLL="cjsre_tmb")
			return(f)
			
		}	
		if(tmbfct=="f3")
		{
			f = MakeADFun(data=list(m=model_data$imat$nocc,ch=model_data$imat$chmat,freq=model_data$imat$freq,frst=model_data$imat$first,
							lst=model_data$imat$last,loc=model_data$imat$loc,tint=model_data$time.intervals,
							phi_fixedDM=phidm,phi_nre=phi_nre,phi_krand=phi_krand,phi_randDM=phi_randDM,
							phi_randIndex=phi_randIndex,phi_counts=phi_counts,phi_idIndex=phi_idIndex,
							p_fixedDM=pdm,p_nre=p_nre,p_krand=p_krand,p_randDM=p_randDM,
							p_randIndex=p_randIndex,p_counts=p_counts,p_idIndex=p_idIndex,getreals=as.integer(getreals),
							prior=as.numeric(prior),mu_phi_prior=mu_phi_prior,sigma_phi_prior=sigma_phi_prior,
							mu_p_prior=mu_p_prior,sigma_p_prior=sigma_p_prior,random_mu_phi_prior=random_mu_phi_prior,
							random_sigma_phi_prior=random_sigma_phi_prior,random_mu_p_prior=mu_p_prior,random_sigma_p_prior=sigma_p_prior),
					parameters=list(phi_beta=initial$Phi,p_beta=initial$p,
							log_sigma_phi=initial$log_sigma_phi,log_sigma_p=rep(-1,npsigma),u_phi=rep(0,phi_nre),u_p=rep(0,p_nre)),
					map=list(phi_beta=factor(rep(NA,length(initial$Phi))),p_beta=factor(rep(NA,length(initial$p))),
															log_sigma_phi=factor(rep(NA,nphisigma)),
															log_sigma_p=factor(rep(NA,npsigma)))
					,DLL="cjsre_tmb")
			return(f)
			
		}	
		
}




