#' Fitting function for Multivariate Multistate CJS with uncertainty models
#' 
#' A function for computing MLEs for MVMS models following Johnson et al 
#' via ADMB. It works very much like mscjs but with more parameters.
#'  
#' It is easiest to call this function through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, Phi.dm,p.dm,Psi.dm,Phi.fixed,p.fixed,Psi.fixed and time.intervals. It is used to save values
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
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param sup supplemental index values for constructing mvms model
#' @param ... not currently used
#' @export
#' @import R2admb
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
#' @references Johnson, D. S., J. L. Laake, S. R. Melin, and DeLong, R.L. 2015. Multivariate State Hidden Markov Models for Mark-Recapture Data. 31:233–244.
mvmscjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		re=FALSE,compile=FALSE,extra.args="",clean=TRUE,sup,...)
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
	chmat=x$ehmat	
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
#	chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=length(ch))
#	for(nlabel in 1:length(strata.labels))
#		chmat=t(apply(chmat,1,sub,pattern=strata.labels[nlabel],replacement=nlabel))
#  Use specified initial values or create if null
	if(is.null(initial))
	{
		par=list(Psi=rep(0,ncol(dml$Psi$fe)),
				p=rep(0,ncol(dml$p$fe)),
				Phi=rep(0,ncol(dml$Phi$fe)))
		if(ncol(dml$pi$fe)>0) par$pi=rep(0,ncol(dml$pi$fe))
		if(ncol(dml$delta$fe)>0) par$delta=rep(0,ncol(dml$delta$fe))
	}
	else
		par=set.initial(names(dml),dml,initial)$par
#  Create list of model data for optimization
	model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe,Psi.dm=dml$Psi$fe,delta.dm=dml$delta$fe,pi.dm=dml$pi$fe,imat=imat,Phi.fixed=parameters$Phi$fixed,
			p.fixed=parameters$p$fixed,Psi.fixed=parameters$Psi$fixed,delta.fixed=parameters$delta$fixed,
			pi.fixed=parameters$pi$fixed,time.intervals=time.intervals)
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
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#  Scale the design matrices and parameters with either input scale or computed scale
	scale=1
	scale=set.scale(names(dml),model_data,scale)
	model_data=scale.dm(model_data,scale)
	# setup tpl to be multistate.tpl 
	if(!re)
		tpl="mvms"
	else
		stop("random effect portion not completed for this model")
	# setup admb exe and cleanup old files and previous tpl; checks for exe 
	# in package directory and if found uses it otherwise compiles tpl file
	setup_admb(tpl,compile,clean,re=FALSE)
	# create .dat file to write its contents 
	con=file(paste(tpl,".dat",sep=""),open="wt")
	# Number of observations
	n=length(model_data$imat$freq)
	write(n,con,append=FALSE)
	# Number of occasions
	nocc=model_data$imat$nocc
	write(nocc,con,append=TRUE)
	# Number of states
	nS=length(strata.labels)
	write(nS,con,append=TRUE)
	# Number of possible observations
    nobs=length(sup$obslevels)
	write(nobs,con,append=TRUE)
	# Number of delta data records per id-occasion
	write(sup$np,con,append=TRUE)
	
	# capture history matrix
	write(t(chmat),con,ncolumns=nocc,append=TRUE)
	# first occasions seen 
	write(model_data$imat$first,con,ncolumns=n,append=TRUE)
	# frequency of capture history 
	if(!re)
	{
		write(model_data$imat$freq,con,ncolumns=n,append=TRUE)
	} else
	{
		if(any(model_data$imat$freq!=1))stop("\n cannot use random effects with frequency >1")
	}
	write(t(model_data$time.intervals),con,ncolumns=nocc-1,append=TRUE)
	# Phi design matrix
	phidm=model_data$Phi.dm
	phifix=rep(-1,nrow(phidm))
	if(!is.null(ddl$Phi$fix))
		phifix[!is.na(ddl$Phi$fix)]=ddl$Phi$fix[!is.na(ddl$Phi$fix)]
	slist=simplify_indices(cbind(phidm,phifix))
	write(ncol(phidm),con,append=TRUE)
	write(length(slist$set),con,append=TRUE)
	write(t(phidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(phidm),append=TRUE)
	write(phifix[slist$set],con,append=TRUE)
	write(slist$indices[ddl$Phi.indices],con,append=TRUE)
	# p design matrix
	pdm=model_data$p.dm
	pfix=rep(-1,nrow(pdm))
	if(!is.null(ddl$p$fix))
		pfix[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
	slist=simplify_indices(cbind(pdm,pfix))
	write(ncol(pdm),con,append=TRUE)
	write(length(slist$set),con,append=TRUE)
	write(t(pdm[slist$set,,drop=FALSE]),con,ncolumns=ncol(pdm),append=TRUE)
	write(pfix[slist$set],con,append=TRUE)
	write(slist$indices[ddl$p.indices],con,append=TRUE)
	
    # write out indices for completing dmat
    write(nrow(sup$indices_forp),con,append=TRUE)
    write(t(sup$indices_forp),con,append=TRUE)

	# Psi design matrix
	# zero out subtracted stratum and remove any unneeded columns
	model_data$Psi.dm[!is.na(ddl$Psi$fix),]=0
	select=vector("logical",length=ncol(model_data$Psi.dm))
	for (i in 1:ncol(model_data$Psi.dm))
		select[i]=any(model_data$Psi.dm[,i]!=0)
	model_data$Psi.dm=model_data$Psi.dm[,select,drop=FALSE]
	psidm=model_data$Psi.dm
	psifix=rep(-1,nrow(psidm))
	if(!is.null(ddl$Psi$fix))
		psifix[!is.na(ddl$Psi$fix)]=ddl$Psi$fix[!is.na(ddl$Psi$fix)]
	slist=simplify_indices(cbind(psidm,psifix))
	write(ncol(psidm),con,append=TRUE)
	write(length(slist$set),con,append=TRUE)
	write(t(psidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(psidm),append=TRUE)
	write(psifix[slist$set],con,append=TRUE)
	write(slist$indices[ddl$Psi.indices],con,append=TRUE)
	# delta design matrix
    deltadm=model_data$delta.dm
    deltafix=rep(-1,nrow(deltadm))
    if(!is.null(ddl$delta$fix))
	   deltafix[!is.na(ddl$delta$fix)]=ddl$delta$fix[!is.na(ddl$delta$fix)]
    slist=simplify_indices(cbind(deltadm,deltafix))
    write(ncol(deltadm),con,append=TRUE)
    write(length(slist$set),con,append=TRUE)
    write(t(deltadm[slist$set,,drop=FALSE]),con,ncolumns=ncol(deltadm),append=TRUE)
    write(deltafix[slist$set],con,append=TRUE)
    write(slist$indices[ddl$delta.indices],con,append=TRUE)
    # pi design matrix
    pidm=model_data$pi.dm
    pifix=rep(-1,nrow(pidm))
    if(!is.null(ddl$pi$fix))
	   pifix[!is.na(ddl$pi$fix)]=ddl$pi$fix[!is.na(ddl$pi$fix)]
    slist=simplify_indices(cbind(pidm,pifix))
    write(ncol(pidm),con,append=TRUE)
    write(length(slist$set),con,append=TRUE)
	if(ncol(pidm)>0)
		write(t(pidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(pidm),append=TRUE)
    write(pifix[slist$set],con,append=TRUE)
    write(slist$indices[ddl$pi.indices],con,append=TRUE)

	if(any(is.na(ddl$pi$fix)))
	   write(0,con,append=TRUE)
	else
	   write(1,con,append=TRUE)
   if(!debug)
	   write(0,con,append=TRUE)
   else
	   write(1,con,append=TRUE)
   
	close(con)
#   write out initial values for betas
	con=file(paste(tpl,".pin",sep=""),open="wt")
	write(par$Phi,con,ncolumns=length(par$Phi),append=FALSE)
	write(par$p,con,ncolumns=length(par$p),append=TRUE)
	write(par$Psi,con,ncolumns=length(par$Psi),append=TRUE)
	if(ncol(dml$pi$fe)>0) 
		write(par$pi,con,ncolumns=length(par$pi),append=TRUE)
	if(ncol(dml$delta$fe)>0) 
	    write(par$delta,con,ncolumns=length(par$delta),append=TRUE)
	close(con)   
	if(hessian)
		xx=run_admb(tpl,extra.args=extra.args)
	else
		xx=run_admb(tpl,extra.args=paste(extra.args,"-nohess"))
	convergence=attr(xx,"status")
	if(is.null(convergence))convergence=0
	res=read_admb(tpl)
	beta=list(unscale.par(c(res$coeflist$phibeta,res$coeflist$pbeta,res$coeflist$dbeta,res$coeflist$psibeta,res$coeflist$pi),scale))
	parnames=names(unlist(beta))
	fixed.npar=length(unlist(beta))
	if(!is.null(res$hes))
	{
		beta.vcv=solvecov(res$hes)$inv
		rownames(res$hes)=parnames
		colnames(res$hes)=rownames(res$hes)
		if(all(diag(beta.vcv>0))) 
			res$cor=beta.vcv/outer(sqrt(diag(beta.vcv)),sqrt(diag(beta.vcv)))
	}  else
		beta.vcv=res$vcov
	rownames(beta.vcv)=parnames
	colnames(beta.vcv)=rownames(beta.vcv)
	rownames(res$cor)=rownames(beta.vcv)
	colnames(res$cor)=rownames(beta.vcv)
	res$vcov=NULL
	optim.details=c(fn=res$fn,maxgrad=res$maxgrad,eratio=res$eratio)
	options=list(extra.args=extra.args)
	res$cor=NULL
	res$maxgrad=NULL
	results=c(beta=beta,neg2lnl=-2*res$loglik,AIC=-2*res$loglik+2*res$npar,convergence=convergence)
	results$optim.details=optim.details
	results$options=options
	results$coeflist=res$coeflist
	results$npar=list(npar=res$npar,npar_sdrpt=res$npar_sdrpt,npar_total=res$npar_total)
	results$beta.vcv=beta.vcv
	res=results
#  Restore non-accumulated, non-scaled dm's etc
	res$model_data=model_data.save
#  Assign S3 class values and return
	class(res)=c("crm","admb","mle","mscjs")
	return(res)
}



