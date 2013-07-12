#' Hidden Markov Model likelihood functions
#' 
#' Functions HMMLikelihood and loglikelihood which compute the log-likelihood for
#' a single capture-recapture sequence and for a set of sequences respectively. The
#' function loglikelihood is called from optimizer and it in turn calls HMMLikelihood for
#' each sequence in the x matrix. These functions are currently limited to a likelihood that
#' is conditional on first capture/release in a known state.
#'  
#' @param id sequential id for the observed sequence
#' @param x single observed sequence (capture history) in HMMLikelihood and matrix of observed sequences (row:id; column:occasion/time) in loglikelihood
#' @param start vector of 2 elements: 1) first state, 2) first occasion or matrix for all id's
#' @param m number of states
#' @param T number of occasions; sequence length
#' @param dmat observation probability matrices
#' @param gamma transition matrices
#' @param par vector of parameter values for log-likelihood evaluation
#' @param type vector of parameter names used to split par vector into types
#' @param freq vector of history frequencies or 1 
#' @param fct_dmat function to create D from parameters
#' @param fct_gamma function to create gamma - transition matrix
#' @param ddl design data list of parameters for each id
#' @param parameters formulas for each parameter type
#' @param debug if TRUE, print out par values and -log-likelihood
#' @usage HMMLikelihood(id,x,start,m,T,dmat,gamma)
#'        loglikelihood(par,type,x,start,m,T,freq=1,fct_dmat,fct_gamma,ddl,parameters,debug=FALSE)
#' @aliases HMMLikelihood loglikelihood
#' @return HMMLikelihood returns log-likelihood for a single sequence and
#' loglikelihood returns the negative log-likelihood for all of the data.
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
HMMLikelihood=function(id,x,start,m,T,dmat,gamma)
{  
	# Arguments:
	# id: sequential id for the observed sequence
	# x: observed sequence (capture (encounter) history)
	# start: vector of 2 elements: 1) first state, 2) first occasion
	# m: number of states
	# T: number of occasions; sequence length
	# parlist: list of parameter vectors split by type (eg Phi, p in CJS)
	# fct_dmat: function to create D from parameters
	# fct_gamma: function to create gamma - transition matrix
	# ddl: design data list of parameters for each id
	# model: formulas for each parameter type
	# Other variables:
	# l: log likelihood value
	# dmat: array of D matrices; one per id and occasion 1 to T-1
	# gamma: array of gamma matrices; one per id and occasion 1 to T-1
	# phi: alpha/sum(alpha) sequence as defined in Zucchini/MacDonald
	# v: temp variable to hold phi calculations
	# u: sum(v)
	l=0
	if(start[2]<T) # compute if first encounter<T 
	{
		# Call specified functions to create dmat and gamma
		# Assign initial prob state vector; phi(x_t) 
		phi=rep(0,m) 
		phi[start[1]]=1
		# Loop over occasions for this encounter history (x)
		for(t in (start[2]+1):T)
		{
			# Compute likelihood contribution for this occasion
			v=phi%*%gamma[t-1,,]%*%diag(dmat[t-1,x[t],])  
			u=sum(v)
			l=l+log(u)
			# Compute updated state vector
			phi=v/u
		}
	}
	return(l)
}     
loglikelihood=function(par,type,x,start,m,T,freq=1,fct_dmat,fct_gamma,ddl,parameters,debug=FALSE)
{
	# Arguments:
	# par: vector of parameter values for log-likelihood evaluation
	# type: vector of parameter names used to split par vector into types
	# x: matrix of observed sequences (row:id; column:occasion/time)
	# start: matrix of initial values; row:id,
	#          2 columns:1) first state, 2) first occasion
	# m: number of states
	# T: number of occasions; sequence length
	# freq: vector of history frequencies or 1 
	# fct_dmat: function to create D from parameters
	# fct_gamma: function to create gamma - transition matrix
	# ddl: design data list of parameters for each id
	# model: formulas for each parameter type
	# Other variables:
	# parlist: list of parameter vectors split by type (eg Phi, p in CJS)
	# gamma: array of transition matrices - one for each id, time
	# dmat: array of observation probability matrices - one for each id, time
	#
	# Create list of parameter vectors from single input parameter vector
	# They are split based on the ordering in the model list; then compute real parameter values
	parlist=split(par,type)
	pars=list()
	for(parname in names(parameters))
		pars[[parname]]=laply(split(reals(parname,ddl=ddl[[parname]],parameters=parameters,parlist=parlist),ddl[[parname]]$id),function(x) x)
    # compute arrays of observation and transition matrices using parameter values
	dmat=fct_dmat(pars,m,T)
	gamma=fct_gamma(pars,m,T)
	# loop over each encounter history in sapply and 
	# create log-likelihood vector - an element for each x
	# sum is total log-likelihood across individuals 
	# return negative log-likelihood
	neglnl=-sum(freq*sapply(1:nrow(x),function(id) 
				 HMMLikelihood(id,x[id,],start[id,],m,T,dmat=dmat[id,,,],gamma=gamma[id,,,])))
	if(debug) cat("\npar = ",par," -lnl= ",neglnl)
	return(neglnl)
}
