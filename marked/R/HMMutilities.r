#' Compute HMM matrices
#' 
#' Computes gamma, dmat and delta (initial) matrices(arrays) and returns them in a list.
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @return list with gamma, dmat and delta arrays
#' @author Jeff Laake
#' @export compute_matrices
#' @keywords utility
compute_matrices=function(object,ddl=NULL)
{
	if(!substr(object$model,1,3)=="HMM")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(is.null(ddl))
		ddl=make.design.data(object$data,object$design.parameters)
	ddl=set.fixed(ddl,object$model.parameters)
	dml=create.dml(ddl,model.parameters=object$model.parameters,design.parameters=object$design.parameters,chunk_size=object$results$options$chunk_size)
	mat=HMMLikelihood(par=object$results$beta,xx=object$data$ehmat,mx=object$data$m,T=object$data$nocc,xstart=object$data$start,freq=object$data$freq,
			fct_dmat=object$data$fct_dmat,fct_gamma=object$data$fct_gamma,fct_delta=object$data$fct_delta,ddl=ddl,dml=dml,parameters=object$model.parameters,
			return.mat=TRUE)
	return(mat)
}
#' Computes backward probabilities 
#' 
#' Computes backward probability sequence for a set of capture histories
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @author Jeff Laake
#' @return array of backward probabilities (one for each id, state, occasion)
#' @export backward_prob
#' @keywords utility
backward_prob=function(object,ddl=NULL)
{  	
	if(!substr(object$model,1,3)=="HMM")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(!is.null(object$results$mat))
	{
		dmat=object$results$mat$dmat
		gamma=object$results$mat$gamma
	}else
	{
		matlist=compute_matrices(object=object,ddl=ddl)
		dmat=matlist$dmat
		gamma=matlist$gamma
	}
	x=object$data$ehmat
	T=object$data$nocc
	first=object$data$start[,2]
	m=object$data$m
	beta=array(NA,dim=c(nrow(x),ncol(x),m))
	# Loop over capture histories
	for(i in 1:nrow(x))
	{
		occ=T
		beta[i,occ,]=rep(1,m)
		# Loop over occasions for this encounter history (x)
		for(t in T:(first[i]+1))
		{
			occ=occ-1
			# Compute backward probability for this occasion
			beta[i,occ,]=gamma[i,t-1,,]%*%diag(dmat[i,t,x[i,t],])%*%beta[i,occ+1,]  
		}
	}
	return(beta)
} 
#' Local decoding of HMM 
#' 
#' Computes state predictions one at a time for each occasion
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @param state.names names for states used to label output; if NULL uses strata.labels + Dead state
#' @author Jeff Laake
#' @return matrix of state predictions
#' @export local_decode
#' @keywords utility
local_decode=function(object,ddl=NULL,state.names=NULL)
{  	
	if(!substr(object$model,1,3)=="HMM")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(is.null(state.names))
		if(!is.null(object$data$strata.labels))
	         state.names=c(object$data$strata.labels,"Dead")
	    else
			 state.names=c("Alive","Dead")
	result=loglikelihood(object)
	result$beta=backward_prob(object)
	stateprob=result$alpha*result$beta/exp(result$lnl)
	states=apply(stateprob,c(1,2),function(x){ if(any(is.na(x))) return(NA) else return(state.names[which.max(x)])})
	return(states)
} 

