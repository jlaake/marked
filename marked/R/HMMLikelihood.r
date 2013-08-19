#' Hidden Markov Model likelihood functions
#' 
#' Functions HMMLikelihood and loglikelihood which compute the log-likelihood for
#' a single capture-recapture sequence and for a set of sequences respectively. The
#' function loglikelihood is called from optimizer and it in turn calls HMMLikelihood for
#' each sequence in the x matrix. 
#'  
#' @param x single observed sequence (capture history) in HMMLikelihood and matrix of observed sequences (row:id; column:occasion/time) in loglikelihood
#' @param first occasion to initiate likelihood calculation for sequence 
#' @param m number of states
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
#' @param parname string name of parameter (eg "Phi")
#' @param parlist list of parameter strings used to split par vector
#' @param start for each ch, the first non-zero x value and the occasion of the first non-zero value
#' @param compute if TRUE, computes reals; otherwise returns column dimension of design matrix for the parameter
#' @usage HMMLikelihood(x,first,m,T,dmat,gamma,delta)
#'        loglikelihood(par,type,x,start,m,T,freq=1,fct_dmat,fct_gamma,fct_delta,ddl,dml,parameters,debug=FALSE)
#'        reals(parname,ddl,dml,parameters,parlist=NULL,compute=TRUE)
#' @aliases HMMLikelihood loglikelihood reals
#' @return HMMLikelihood returns log-likelihood for a single sequence and
#' loglikelihood returns the negative log-likelihood for all of the data. reals
#' returns either the column dimension of design matrix for parameter or the real parameter vector
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
HMMLikelihood=function(x,first,m,T,dmat,gamma,delta)
{  
	# Arguments:
	# x: observed sequence (capture (encounter) history)
	# first: occasion to start sequence
	# m: number of states
	# T: number of occasions; sequence length
	# dmat: array of occasion specific observation probabilty matrices
	# gamma: array of occasion specific transition matrices
	# delta: initial state distribution
	# Other variables:
	# lnl: log likelihood value
	# phi: alpha/sum(alpha) sequence as defined in Zucchini/MacDonald
	# v: temp variable to hold phi calculations
	# u: sum(v)
	# Assign prob state vector for initial observation: delta*p(x_first)
 	v=delta%*%diag(dmat[first,x[first],]) 
	# Compute log-likelihood contribution for first observation; for
	# models that condition on first observation u=1,lnl=0
	u=sum(v)
	phi=v/u
	lnl=log(u)
	# Loop over occasions for this encounter history (x)
	for(t in (first+1):T)
	{
		# Compute likelihood contribution for this occasion
		v=phi%*%gamma[t-1,,]%*%diag(dmat[t,x[t],])  
		u=sum(v)
		lnl=lnl+log(u)
		# Compute updated state vector
		phi=v/u
	}
	return(lnl)
}     
loglikelihood=function(par,type,x,start,m,T,freq=1,fct_dmat,fct_gamma,
		fct_delta,ddl,dml,parameters,debug=FALSE)
{
	# Arguments:
	# par: vector of parameter values for log-likelihood evaluation
	# type: vector of parameter names used to split par vector into types
	# x: matrix of observed sequences (row:id; column:occasion/time)
	# start: matrix with a row for each id and 2 columns 
	#        1) first observed state, 2) first occasion observed
	# m: number of states
	# T: number of occasions; sequence length
	# freq: vector of history frequencies or 1 
	# fct_dmat: function to create D from parameters
	# fct_gamma: function to create gamma - transition matrix
	# fct_delta: function to create initial state probability distribution matrix
	# ddl: design data list of parameters for each id
	# model: formulas for each parameter type
	# Other variables:
	# parlist: list of parameter vectors split by type (eg Phi, p in CJS)
	# gamma: array of transition matrices - one for each id, time
	# dmat: array of observation probability matrices - one for each id, time
	#
	# Create list of parameter matrices from single input parameter vector
	# Firs split parameter vector by prameter type (type) 
	parlist=split(par,type)
	pars=list()
	# For each parameter type call function reals to compute vector
	# of real parameter values; then use laply and split to create
	# a matrix of parameter values with a row for each id and column for
	# each occasion.
    for(parname in names(parameters))
    {
        R=reals(parname,ddl=ddl[[parname]],dml=dml[[parname]],parameters=parameters,parlist=parlist)
        pars[[parname]]=laply(split(R,ddl[[parname]]$id),function(x) x)
    }
	# compute 4-d arrays of id- and occasion-specific 
	#observation and transition matrices using parameter values
	dmat=fct_dmat(pars,m,T)
	gamma=fct_gamma(pars,m,T)
	# compute matrix of initial state distribution for each id
	delta=fct_delta(pars,m,T,start)
#	browser()
#	xx=hmm.lnl(x,start,m,T,dmat,gamma,delta)
	# loop over each encounter history in sapply and 
	# create log-likelihood vector - an element for each x
	# sum is total log-likelihood across individuals 
	# return negative log-likelihood
	neglnl=-sum(freq*sapply(1:nrow(x),function(id) 
						HMMLikelihood(x[id,],start[id,2],m,T,
								dmat=dmat[id,,,],gamma=gamma[id,,,],
								delta=delta[id,])))
	if(debug){
		cat("\npar \n")
		print(split(par,type))
		for(parname in names(parameters))
		{
			cat("\n",parname,"\n")
			print(pars[[parname]][1,])
		}
		cat(" -lnl= ",neglnl)
		ps=delta[1,]
		for(i in 1:(T-1))
		{
			ps=ps%*%gamma[1,i,,]
			cat("\ni = ",i+1," ps= ", ps)
		}
		
	}
	return(neglnl)
}
# Computes real estimates for HMM models using inverse of link from design data (ddl) and model for a particular parameter type (parname) or
# returns the number of columns in the design matrix (compute=FALSE); handles fixed parameters assigned by non-NA value in field named 
# fix in the ddl dataframe.
reals=function(parname,ddl,dml,parameters,parlist=NULL,compute=TRUE)
{
	dm=dml$fe
	#dm=model.matrix(parameters[[parname]]$formula,ddl)
	# if some reals are fixed, assign 0 to rows of dm and then
	# remove any columns (parameters) that are all 0.
	if(!is.null(ddl$fix)&&any(!is.na(ddl$fix)))
	{
		dm[!is.na(ddl$fix),]=0
		dm=dm[,apply(dm,2,function(x) any(x!=0)),drop=FALSE]
	}
	# if not computing reals, return the names of columns in dm
	if(!compute)return(colnames(dm))
	# Currently for log or logit link, return the inverse values
	values=switch(parameters[[parname]]$link,
	     log=exp(as.vector(dm%*%parlist[[parname]])),
		 logit=plogis(as.vector(dm%*%parlist[[parname]])),
		 identity=as.vector(dm%*%parlist[[parname]]))
    if(!is.null(ddl$time.interval))values=values^ddl$time.interval
	# if some reals are fixed, set reals to their fixed values 
	if(!is.null(ddl$fix))
		values[!is.na(ddl$fix)]=ddl$fix[!is.na(ddl$fix)]
	# return vector of reals
	return(values)
}

hmm.lnl=function(x,start,m,T,dmat,gamma,delta)
{
	browser()
	lnl=.Fortran("hmmlike",as.integer(x),as.integer(nrow(x)),as.integer(m),as.integer(T),
			as.integer(nrow(dmat[1,1,,])),as.integer(start),as.double(dmat),
			as.double(gamma),as.double(delta),lnl=double(lnl),P=double(P),PACKAGE="marked")
    lnl
}
		