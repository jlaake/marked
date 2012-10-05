#' Set initial values 
#' 
#' Sets initial values specified in a list.
#' 
#' @param pars character vector of parameter names
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param initial list of vectors for parameter initial values
#' @return List of initial values for each parameter in the model
#' @author Jeff Laake <jeff.laake@@noaa.gov>
set.initial=function(pars,dml,initial)
{
	if(class(initial)[1]=="crm")
		if(class(initial)[2]=="mcmc")
			initial=lapply(initial$results$beta,function(x){z=x$mean;names(z)=rownames(x);z})	    
	    else
			initial=initial$results$beta
	par=vector("list",length(pars))
	names(par)=pars
	for(parx in pars)
	{
		init=initial[[parx]]
		if(is.null(init))init=0
		if(length(init)==1 &is.null(names(init)))
			par[[parx]]=c(init,rep(0,ncol(dml[[parx]])-1))
		else
		{
			if(is.null(names(init)))
			{
				if(length(init)!=ncol(dml[[parx]]))
					stop(paste("For",parx,",length of initial vector does not match number of parameters."))
				else
					par[[parx]]=initial
			} else
			{
				beta.names=colnames(dml[[parx]])
				par[[parx]]=rep(0,length(beta.names))
				par[[parx]][beta.names%in%names(init)]=init[which(names(init)%in%beta.names)]
			}
		}
	}
	return(par)
}