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
	mat=HMMLikelihood(par=object$results$beta,xx=object$data$ehmat,mx=data$m,T=object$data$nocc,xstart=object$data$start,freq=object$data$freq,
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
		matlist=compute_matrices(x=x,ddl=ddl)
		dmat=matlist$dmat
		gamma=matlst$gamma
	}
	x=object$data$ehmat
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
			beta[i,occ,]=gamma[i,t-1,,]%*%diag(dmat[i,t,x[t],])%*%beta[i,occ+1,]  
		}
	}
	return(beta)
} 

