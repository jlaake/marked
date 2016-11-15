#' Hidden Markov Model likelihood functions
#' 
#' Function HMMLikelihood computes the log-likelihood via hmm.lnl which is a wrapper for the
#' FORTRAN code hmm_like.f.  The function HMMlikelihood is called from optimizer and it in turn calls hmm.lnl after
#' setting up parameters.
#'  
#' For an R version of the HMMLikelihood and related code see \code{\link{R_HMMLikelihood}}
#'  
#' @param xx  matrix of observed sequences (row:id; column:occasion/time); xx used instead of x to avoid conflict in optimx
#' @param x  same as xx but for call to hmm.lnl
#' @param mx number of states; mx used instead of m to avoid conflict in optimx
#' @param m  same as mx but for call to hmm.lnl
#' @param T number of occasions; sequence length
#' @param dmat observation probability matrices
#' @param gamma transition matrices
#' @param delta initial distribution
#' @param par vector of parameter values for log-likelihood evaluation
#' @param type vector of parameter names used to split par vector into types
#' @param freq vector of history frequencies or 1 
#' @param fct_dmat function to create D from parameters
#' @param fct_gamma function to create gamma - transition matrix
#' @param fct_delta function to create initial state distribution
#' @param ddl design data list of parameters for each id
#' @param dml list of design matrices; one entry for each parameter; each entry contains fe and re for fixed and random effects
#' @param parameters formulas for each parameter type
#' @param debug if TRUE, print out par values and -log-likelihood
#' @param parlist list of parameter strings used to split par vector
#' @param xstart for each ch, the first non-zero x value and the occasion of the first non-zero value; ; xstart used instead of start to avoid conflict in optimx
#' @param start same as xstart but for hmm.lnl
#' @param return.mat If TRUE, returns list of transition, observation and delta arrays.
#' @param sup  list of supplemental information that may be needed by the function but only needs to be computed once; currently only used for MVMS models for dmat
#' @param check if TRUE, checks validity of gamma, dmat and delta to look for any errors
#' @param indices specific indices for computation unless NULL
#' @usage HMMLikelihood(par,type,xx,xstart,mx,T,freq=1,fct_dmat,fct_gamma,fct_delta,ddl,
#'                          dml,parameters,debug=FALSE,return.mat=FALSE,sup=NULL,check=FALSE)
#'        reals(ddl,dml,parameters,parlist,indices=NULL)
#'        hmm.lnl(x,start,m,T,dmat,gamma,delta,freq,debug)
#' @aliases HMMLikelihood reals hmm.lnl
#' @return HMMLikelihood returns log-likelihood for a single sequence and
#' hmm.lnl returns the negative log-likelihood value for each capture history. reals
#' returns either the column dimension of design matrix for parameter or the real parameter vector
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @seealso R_HMMLikelihood
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
HMMLikelihood=function(par,type=NULL,xx,xstart,mx,T,freq=1,fct_dmat,fct_gamma,
		fct_delta,ddl,dml,parameters,debug=FALSE,return.mat=FALSE,sup=NULL,check=FALSE)
{
	m=mx
	# Create list of parameter matrices from single input parameter vector
	# First split parameter vector by prameter type (type) 
	# changed to mx to avoid problems with optimx
	if(is.null(type))
		parlist=par
	else
		parlist=split(par,type)
	# For each parameter type call function reals to compute vector
	# of real parameter values; then use split to create
	# a matrix of parameter values with a row for each id and column for
	# each occasion.
	pars=list()
	for(parname in names(parameters))
    {
		indices=ddl[[paste(paste(parname,"indices",sep="."))]]
        R=reals(ddl=ddl[[parname]],dml=dml[[parname]],parameters=parameters[[parname]],parlist=parlist[[parname]],indices=indices)
		R[is.infinite(R)]=1e199*sign(Inf)
        if(is.null(indices))
		  pars[[parname]]=do.call("rbind",split(R,ddl[[parname]]$id))
	  else
		  pars[[parname]]=do.call("rbind",split(R,ddl[[paste(parname,"id",sep=".")]]))
  }
	# compute 4-d arrays of id- and occasion-specific 
	# observation probability matrices 
	if(is.null(sup))
       dmat=fct_dmat(pars,m,F=xstart[,2],T)
   else
   {
	   # unkinit set to 1 unless all unknown or all known at initial release;
	   # when unkinit=1, delta is applied in dmat
	   sup$unkinit=1-as.numeric(all(is.na(xstart[,1])) | all(!is.na(xstart[,1])))
	   dmat=fct_dmat(pars,m,F=xstart[,2],T,sup)
   }
   if(check)
   {
	   bad=FALSE
	   for(i in 1:dim(dmat)[1])
	   for(j in (xstart[i,2]+1):dim(dmat)[2])
	   for(k in 1:dim(dmat)[4])
		   if(is.nan(sum(dmat[i,j,,k])) || (sum(dmat[i,j,,k])< -1e-8 | sum(dmat[i,j,,k])> 1.0000001))
		   {
			  cat("\nFor id = ",i," occasion = ", j, " and strata = ",k," value sum of dmat is",sum(dmat[i,j,,k]))
			  bad=TRUE
		  }
	   if(bad) stop("\nCheck model and design data fix values for parameters affecting dmat -the observation probability matrix\n")   
   }
   # transition matrices using parameter values
	gamma=fct_gamma(pars,m,F=xstart[,2],T)
	if(check)
	{
		bad=FALSE
		for(i in 1:dim(gamma)[1])
			for(j in 1:dim(gamma)[2])
				for(k in 1:dim(gamma)[3])
					if(is.nan(sum(gamma[i,j,k,])) || (sum(gamma[i,j,k,])< -1e-8 | sum(gamma[i,j,k,])> 1.0000001))
					{
						cat("\nFor id = ",i," occasion = ", j, " and strata = ",k," value sum of gamma is",sum(gamma[i,j,k,]))
						bad=TRUE
					}
		if(bad)stop("\nCheck model and design data fix values for parameters affecting gamma -the transition probability matrix\n")
	}
	
	# compute matrix of initial state distribution for each id
	delta=fct_delta(pars,m,F=xstart[,2],T,xstart)
	if(check)
	{
		zz=rowSums(delta)
	    if(any(is.nan(zz)) || any((zz< (1-1e-8) | zz> (1+1e-8) )))
	    {
		   cat("\nProblem with initial state for observations",paste((1:nrow(delta))[is.nan(zz) || (zz< (1-1e-8) | zz> (1+1e-8))],sep=","))
		   stop("\nCheck model and design data fix values for parameters affecting delta - initial state probability")
	    }	
	}
	if(return.mat)	
		return(list(dmat=dmat,gamma=gamma,delta=delta))
	if(is.list(m)) m=m$ns*m$na+1
	if(debug){
		cat("\npar \n")
		print(parlist)
	}
	neglnl=hmm.lnl(xx,xstart[,2],m,T,dmat,gamma,delta,rowSums(freq),debug)
	if(debug)cat("\n -lnl= ",neglnl)
	return(neglnl)
}
reals=function(ddl,dml,parameters,parlist,indices=NULL)
{
	# Computes real estimates for HMM models using inverse of 
	# link from design matrix and for a particular parameter 
	# type (parname); handles fixed parameters assigned by 
	# non-NA value in field named fix in the ddl dataframe.
	dm=dml$fe
	if(is.null(dm))return(NULL)
	if(ncol(dm)!=0)
	{	
		# Currently for log,logit or identity link, return the inverse values
		values=switch(parameters$link,
				log=exp(as.vector(dm%*%parlist)),
				logit=plogis(as.vector(dm%*%parlist)),
				identity=as.vector(dm%*%parlist))
		if(!is.null(indices))
			values = values[indices]
		#if(!is.null(dml$indices))values=values[dml$indices]
		if(!is.null(ddl$time.interval))values=values^ddl$time.interval
	}
	else
	{
		values=rep(NA,nrow(dm))
		if(!is.null(indices))
			values = values[indices]
	}
	# if some reals are fixed, set reals to their fixed values 
	if(!is.null(ddl$fix)) 
	{
		if(is.null(indices))
			values[!is.na(ddl$fix)]=ddl$fix[!is.na(ddl$fix)]
		else
		{
			fix=ddl$fix[indices]
			values[!is.na(fix)]=fix[!is.na(fix)]
		}
	}		
	# return vector of reals
	return(values)
}
hmm.lnl=function(x,start,m,T,dmat,gamma,delta,freq,debug)
{
	lnl=.Fortran("hmmlike",as.integer(x),as.integer(nrow(x)),as.integer(m),as.integer(T),
			as.integer(nrow(dmat[1,1,,])),as.integer(start),as.double(freq),as.double(dmat),
			as.double(gamma),as.double(delta),lnl=double(nrow(x)),PACKAGE="marked")$lnl
	if(any(is.nan(lnl) | is.infinite(lnl)))
	{
		if(debug)
		{
			cat("lnl is nan for these rows\n")
		    cat(which(is.nan(lnl)))
			cat("lnl is infinite for these rows\n")
			cat(which(is.infinite(lnl)))
		}
	    lnl[is.nan(lnl)]=-1e5
		lnl[is.infinite(lnl)]=-1e5
	}
    -sum(lnl)
}
