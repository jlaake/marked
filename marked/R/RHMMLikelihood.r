#' Hidden Markov Model Functions
#' 
#' R implementation of HMMs described in processed report. Function HMMLikelihood renamed to R_HMMLikelihood and
#' loglikelihood modified to work with a single x and to return lnl, alpha, phi, v, dmat,  and gamma values. backward_prob returns the backward probabilities.
#' These are not used by the fitting code.
#'  
#' @param x single observed sequence (capture history) 
#' @param id id of capture history for computation
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
#' @param start for each ch, the first non-zero x value and the occasion of the first non-zero value
#' @usage R_HMMLikelihood(x,first,m,T,dmat,gamma,delta)
#'        loglikelihood(par,type,x,id,start,m,T,freq=1,fct_dmat,fct_gamma,fct_delta,
#'                         ddl,dml,parameters)
#'        backward_prob(x,first,m,T,dmat,gamma)
#' @aliases R_HMMLikelihood loglikelihood backward_prob
#' @return both return log-likelihood, alpha and phi vectors and dmat and gamma matrices for a single sequence
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
R_HMMLikelihood=function(x,first,m,T,dmat,gamma,delta)
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
	alpha=matrix(NA,nrow=T-first+1,ncol=m)
	phimat=matrix(NA,nrow=T-first+1,ncol=m)
	vmat=matrix(NA,nrow=T-first+1,ncol=m)
	occ=1
	# Assign prob state vector for initial observation: delta*p(x_first)
	v=delta%*%diag(dmat[first,x[first],]) 
	# Compute log-likelihood contribution for first observation; for
	# models that condition on first observation u=1,lnl=0
	u=sum(v)
	phi=v/u
	alpha[occ,]=v
	vmat[occ,]=v
	phimat[occ,]=phi
	lnl=log(u)
	# Loop over occasions for this encounter history (x)
	for(t in (first+1):T)
	{
		occ=occ+1
		# Compute likelihood contribution for this occasion
		v=phi%*%gamma[t-1,,]%*%diag(dmat[t,x[t],])  
		u=sum(v)
		lnl=lnl+log(u)
		# Compute updated state vector
		phi=v/u
		# Compute alpha vector
		vmat[occ,]=v
		alpha[occ,]=alpha[occ-1,]%*%gamma[t-1,,]%*%diag(dmat[t,x[t],])
		phimat[occ,]=phi
	}
	beta=backward_prob(x,first,m,T,dmat,gamma)
	dmat=apply(dmat,c(2,3,1),function(x)x)
	gamma=apply(gamma,c(2,3,1),function(x)x)
	dimnames(dmat)[3]=list(1:T)
	dimnames(gamma)[3]=list(1:(T-1))
	dimnames(alpha)[1]=list(first:T)
	dimnames(beta)[1]=list(first:T)
	dimnames(phimat)[1]=dimnames(alpha)[1]
	dimnames(vmat)[1]=dimnames(alpha)[1]
	return(list(lnl=lnl,alpha=alpha,beta=beta,phi=phimat,v=vmat,dmat=dmat,gamma=gamma))
} 
loglikelihood=function(par,type,x,id,start,m,T,freq=1,fct_dmat,fct_gamma,
		fct_delta,ddl,dml,parameters)
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
	# First split parameter vector by prameter type (type) 
	parlist=split(par,type)
	pars=list()
	# For each parameter type call function reals to compute vector
	# of real parameter values; then use laply and split to create
	# a matrix of parameter values with a row for each id and column for
	# each occasion.
	for(parname in names(parameters))
	{
		R=reals(ddl=ddl[[parname]],dml=dml[[parname]],parameters=parameters[[parname]],parlist=parlist[[parname]])
		pars[[parname]]=laply(split(R,ddl[[parname]]$id),function(x) x)
	}
	# compute 4-d arrays of id- and occasion-specific 
	#observation and transition matrices using parameter values
	dmat=fct_dmat(pars,m,F=start[,2],T)
	gamma=fct_gamma(pars,m,F=start[,2],T)
	# compute matrix of initial state distribution for each id
	delta=fct_delta(pars,m,F=start[,2],T,start)
	# loop over each encounter history in sapply and 
	# create log-likelihood vector - an element for each x
	# sum is total log-likelihood across individuals 
	# return negative log-likelihood
	object=R_HMMLikelihood(x[id,],start[id,2],m,T,
			dmat=dmat[id,,,],gamma=gamma[id,,,],
			delta=delta[id,])
#	object$beta=backward_prob(x[id,],start[id,2],m,T,
#								dmat=dmat[id,,,],gamma=gamma[id,,,])
	return(object)
}

