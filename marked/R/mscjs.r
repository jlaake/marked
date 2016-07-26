#' Fitting function for Multistate CJS models
#' 
#' A function for computing MLEs for a Multi-state Cormack-Jolly-Seber open
#' population capture-recapture model for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{mscjs} through the function \code{\link{crm}}.
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
#' @references Ford, J. H., M. V. Bravington, and J. Robbins. 2012. Incorporating individual variability into mark-recapture models. Methods in Ecology and Evolution 3:1047-1054.
#' @examples 
#' \donttest{
#' # this example requires admb
#' # The same example is in the RMark package and it is included here to
#' # illustrate the differences in the handling of mlogit parameters between RMark 
#' # and marked.  The MARK software handles parameters like Psi which must sum to 1
#' # by excluding one of the cells that is used as a reference cell and is computed by
#' # subtracting the other cell values from 1 so the total sums to 1.  This is often
#' # handled with an mlogit parameter in which the cell values are exp(beta) and the
#' # reference cell is set to 1 and the values are divided by the sum across the cells
#' # so the resulting values are probabilities that sum to 1. In marked, instead of removing
#' # one of the cells, all are included and the user must select which should be the
#' # reference cell by setting the value fix=1 for that cell and others are NA so they are
#' # estimated. For transition parameters like Psi, the default design data is setup so 
#' # that the probability of remaining in the cell (stratum=tostratum) is the reference cell
#' # and fix set to 1.  Thus, this means 2 changes are needed to the script in RMark.
#' # The first is to remove the statement skagit.ddl$Psi$fix=NA because that over-rides
#' # the default fix values.  The other is to add
#' # skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="B"&
#' #  skagit.ddl$Psi$time==5]=0
#' # to change the value from 1 to 0 which forces movement from B to A in the interval 5 to 6. If
#' # this is not done then Psi B to B=Psi B to A=0.5 because each is 1 and when they are normalized
#' # they are divided by the sum which is 2 (1/2).
#' if(class(try(setup_admb("mscjs")))!="try-error")
#' {
#' data(skagit)
#' skagit.processed=process.data(skagit,model="Mscjs",groups=c("tag"),strata.labels=c("A","B"))
#'skagit.ddl=make.design.data(skagit.processed)
#'#
#'# p
#'#
#'# Can't be seen at 5A or 2B,6B (the latter 2 don't exist)
#'skagit.ddl$p$fix=ifelse((skagit.ddl$p$stratum=="A"&skagit.ddl$p$time==5) | 
#' (skagit.ddl$p$stratum=="B"&skagit.ddl$p$time%in%c(2,6)),0,NA)
#'# Estimated externally from current data to allow estimation of survival at last interval
#'skagit.ddl$p$fix[skagit.ddl$p$tag=="v7"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.687
#'skagit.ddl$p$fix[skagit.ddl$p$tag=="v9"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.975
#'#
#'# Psi
#'#
#'# only 3 possible transitions are A to B at time interval 2 to 3 and 
#'# for time interval 3 to 4 from A to B and from B to A
#'# rest are fixed values
#'############ change for RMark to marked; remove next line
#'#skagit.ddl$Psi$fix=NA
#'# stay in A for intervals 1-2, 4-5 and 5-6
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="A"&
#'  skagit.ddl$Psi$tostratum=="B"&skagit.ddl$Psi$time%in%c(1,4,5)]=0
#'# stay in B for interval 4-5
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"
#'  &skagit.ddl$Psi$time==4]=0
#'# leave B to go to A for interval 5-6
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
#' skagit.ddl$Psi$time==5]=1
#'############ change for RMark to marked; add next line to set B to B to 0 otherwise it has
#'############ been set to 1 by default which would make psi B to B = psi B to A = 0.5
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="B"&
#' skagit.ddl$Psi$time==5]=0
#'# "stay" in B for interval 1-2 and 2-3 because none will be in B
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
#' skagit.ddl$Psi$time%in%1:2]=0
#'# 
#'# S
#'#
#'# None in B, so fixing S to 1
#'skagit.ddl$S$fix=ifelse(skagit.ddl$S$stratum=="B"&skagit.ddl$S$time%in%c(1,2),1,NA)
#'skagit.ddl$S$fix[skagit.ddl$S$stratum=="A"&skagit.ddl$S$time==4]=1
#'# fit model
#'p.timexstratum.tag=list(formula=~time:stratum+tag,remove.intercept=TRUE)
#'Psi.sxtime=list(formula=~-1+stratum:time)
#'S.stratumxtime=list(formula=~-1+stratum:time)
#'#
#'mod1=crm(skagit.processed,skagit.ddl,
#' model.parameters=list(S=S.stratumxtime,p= p.timexstratum.tag,Psi=Psi.sxtime),hessian=TRUE)
#' if(class(mod1)[1]!="try-error") mod1
#'} }
mscjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
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
	# setup tpl to be multistate.tpl 
	if(!re)
		tpl="multistate"
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
	# S design matrix
    phidm=as.matrix(model_data$S.dm)
    phifix=rep(-1,nrow(phidm))
	if(!is.null(ddl$S$fix))
	   phifix[!is.na(ddl$S$fix)]=ddl$S$fix[!is.na(ddl$S$fix)]
    #slist=simplify_indices(cbind(phidm,phifix))
    write(ncol(phidm),con,append=TRUE)
	write(nrow(phidm),con,append=TRUE)
	write(t(phidm),con,ncolumns=ncol(phidm),append=TRUE)
	write(phifix,con,append=TRUE)
	write(ddl$S.indices,con,append=TRUE)
#	write(length(slist$set),con,append=TRUE)
#	write(t(phidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(phidm),append=TRUE)
#	write(phifix[slist$set],con,append=TRUE)
#	write(slist$indices,con,append=TRUE)
	# p design matrix
    # zero out dm if any unobserved stratum; only done to remove unneeded columns
#    if(uS>0)
#	{
#		model_data$p.dm[as.numeric(ddl$p$stratum)>=(length(strata.labels)-uS),]=0
#		select=vector("logical",length=ncol(model_data$p.dm))
#		for (i in 1:ncol(model_data$p.dm))
#			select[i]=any(model_data$p.dm[,i]!=0)
#		model_data$p.dm=model_data$p.dm[,select,drop=FALSE]
#	}
	pdm=as.matrix(model_data$p.dm)
	pfix=rep(-1,nrow(pdm))
	if(!is.null(ddl$p$fix))
		pfix[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
	slist=simplify_indices(cbind(pdm,pfix))
	write(ncol(pdm),con,append=TRUE)
	write(nrow(pdm),con,append=TRUE)
	write(t(pdm),con,ncolumns=ncol(pdm),append=TRUE)
	write(pfix,con,append=TRUE)
	write(ddl$p.indices,con,append=TRUE)
#	write(length(slist$set),con,append=TRUE)
#	write(t(pdm[slist$set,,drop=FALSE]),con,ncolumns=ncol(pdm),append=TRUE)
#	write(pfix[slist$set],con,append=TRUE)
#	write(slist$indices,con,append=TRUE)
	
	# Psi design matrix
	# zero out subtracted stratum and remove any unneeded columns
#    model_data$Psi.dm[!is.na(ddl$Psi$fix),]=0
#	select=vector("logical",length=ncol(model_data$Psi.dm))
#	for (i in 1:ncol(model_data$Psi.dm))
#		select[i]=any(model_data$Psi.dm[,i]!=0)
#	model_data$Psi.dm=model_data$Psi.dm[,select,drop=FALSE]
	psidm=as.matrix(model_data$Psi.dm)
	psifix=rep(-1,nrow(psidm))
	if(!is.null(ddl$Psi$fix))
	psifix[!is.na(ddl$Psi$fix)]=ddl$Psi$fix[!is.na(ddl$Psi$fix)]
    slist=simplify_indices(cbind(psidm,psifix))
	write(ncol(psidm),con,append=TRUE)
	write(nrow(psidm),con,append=TRUE)
	write(t(psidm),con,ncolumns=ncol(psidm),append=TRUE)
	write(psifix,con,append=TRUE)
	write(ddl$Psi.indices,con,append=TRUE)
#    write(length(slist$set),con,append=TRUE)
#   write(t(psidm[slist$set,,drop=FALSE]),con,ncolumns=ncol(psidm),append=TRUE)
#    write(psifix[slist$set],con,append=TRUE)
#    write(slist$indices,con,append=TRUE)
	close(con)
#   write out initial values for betas
	con=file(paste(tpl,".pin",sep=""),open="wt")
	write(par$S,con,ncolumns=length(par$S),append=FALSE)
	write(par$p,con,ncolumns=length(par$p),append=TRUE)
	write(par$Psi,con,ncolumns=length(par$Psi),append=FALSE)
	close(con)   
	if(hessian)
		xx=run_admb(tpl,extra.args=extra.args)
	else
		xx=run_admb(tpl,extra.args=paste(extra.args,"-nohess"))
	convergence=attr(xx,"status")
	if(is.null(convergence))convergence=0
	res=read_admb(tpl)
	beta=list(unscale.par(c(res$coeflist$phibeta,res$coeflist$pbeta,res$coeflist$psibeta),scale))
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


