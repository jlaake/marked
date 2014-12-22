#' Hidden Markov Model likelihood functions
#' 
#' Function HMMLikelihood computes the log-likelihood via hmm.like which is a wrapper for the
#' FORTRAN code hmm_like.f.  The function HMMlikelihood is called from optimizer and it in turn calls HMMLikelihood for
#' each sequence in the x matrix.
#'  
#' For an R version of the HMMLikelihood and related code see \code{\link{RHMMLikelihood}}
#'  
#' @param xx  matrix of observed sequences (row:id; column:occasion/time); xx used instead of x to avoid conflict in optimx
#' @param mx number of states; mx used instead of m to avoid conflict in optimx
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
#' @param return.mat If TRUE, returns list of transition, observation and delta arrays.
#' @usage HMMLikelihood(par,type,xx,xstart,mx,T,freq=1,fct_dmat,fct_gamma,fct_delta,ddl,
#'                          dml,parameters,debug=FALSE,return.mat=FALSE)
#'        reals(ddl,dml,parameters,parlist)
#'        hmm.lnl(x,start,m,T,dmat,gamma,delta,freq)
#' @aliases HMMLikelihood reals hmm.lnl
#' @return HMMLikelihood returns log-likelihood for a single sequence and
#' loglikelihood returns the negative log-likelihood for all of the data. reals
#' returns either the column dimension of design matrix for parameter or the real parameter vector
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @seealso RHMMLikelihood
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
HMMLikelihood=function(par,type=NULL,xx,xstart,mx,T,freq=1,fct_dmat,fct_gamma,
		fct_delta,ddl,dml,parameters,debug=FALSE,return.mat=FALSE)
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
	# of real parameter values; then use laply and split to create
	# a matrix of parameter values with a row for each id and column for
	# each occasion.
	pars=list()
	for(parname in names(parameters))
    {
        R=reals(ddl=ddl[[parname]],dml=dml[[parname]],parameters=parameters[[parname]],parlist=parlist[[parname]])
		R[is.infinite(R)]=1e199*sign(Inf)
        pars[[parname]]=laply(split(R,ddl[[parname]]$id),function(x) x)
    }
	# compute 4-d arrays of id- and occasion-specific 
	# observation probability matrices and 
    dmat=fct_dmat(pars,m,F=xstart[,2],T)
	# transition matrices using parameter values
	gamma=fct_gamma(pars,m,F=xstart[,2],T)
	# compute matrix of initial state distribution for each id
	delta=fct_delta(pars,m,F=xstart[,2],T,xstart)
	if(return.mat)return(list(dmat=dmat,gamma=gamma,delta=delta))
	if(is.list(m)) m=m$ns*m$na+1
	if(debug){
		cat("\npar \n")
		print(parlist)
	}
	neglnl=hmm.lnl(xx,xstart,m,T,dmat,gamma,delta,rowSums(freq),debug)
	if(debug)cat("\n -lnl= ",neglnl)
	return(neglnl)
}
reals=function(ddl,dml,parameters,parlist)
{
	# Computes real estimates for HMM models using inverse of 
	# link from design matrix and for a particular parameter 
	# type (parname); handles fixed parameters assigned by 
	# non-NA value in field named fix in the ddl dataframe.
	dm=dml$fe
	if(ncol(dm)!=0)
	{	
		# Currently for log,logit or identity link, return the inverse values
		values=switch(parameters$link,
				log=exp(as.vector(dm%*%parlist)),
				logit=plogis(as.vector(dm%*%parlist)),
				identity=as.vector(dm%*%parlist))
		if(!is.null(ddl$time.interval))values=values^ddl$time.interval
	}
	else
		values=rep(NA,nrow(dm))
	# if some reals are fixed, set reals to their fixed values 
	if(!is.null(ddl$fix))
		values[!is.na(ddl$fix)]=ddl$fix[!is.na(ddl$fix)]
	# return vector of reals
	return(values)
}
hmm.lnl=function(x,start,m,T,dmat,gamma,delta,freq,debug)
{
	lnl=.Fortran("hmmlike",as.integer(x),as.integer(nrow(x)),as.integer(m),as.integer(T),
			as.integer(nrow(dmat[1,1,,])),as.integer(start),as.double(freq),as.double(dmat),
			as.double(gamma),as.double(delta),lnl=double(nrow(x)),PACKAGE="marked")$lnl
	if(any(is.nan(lnl)))
	{
		if(debug)
		{
			cat("lnl is nan for these rows\n")
		    cat(which(is.nan(lnl)))
		}
	    lnl[is.nan(lnl)]=-1e5
	}
    -sum(lnl)
}
