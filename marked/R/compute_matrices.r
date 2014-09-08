#' Compute HMM matrices
#' 
#' Computes gamma, dmat and intial matrices(arrays) and returns them in a list.
#' 
#' @param x fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @author Jeff Laake
#' @export 
#' @keywords utility
compute_matrices=function(x,ddl)
{
	if(!substr(x$model,1,3)=="HMM")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	m=list(ns=length(x$data$strata.list$states),na=length(x$data$strata.list[[names(x$data$strata.list)[names(x$data$strata.list)!="states"]]]))
	ddl=set.fixed(ddl,x$model.parameters)
	dml=create.dml(ddl,model.parameters=x$model.parameters,design.parameters=x$design.parameters,chunk_size=1e7)
	mat=HMMLikelihood(par=x$results$beta,xx=x$data$ehmat,mx=m,T=x$data$nocc,start=x$data$start,freq=x$data$freq,
			fct_dmat=x$data$fct_dmat,fct_gamma=x$data$fct_gamma,fct_delta=x$data$fct_delta,ddl=ddl,dml=dml,parameters=x$model.parameters,
			return.mat=TRUE)
    return(mat)
}
