#' Scaling functions
#' 
#' Set scale, scale dm and scale/unscale parameters 
#'
#' @usage 	set_scale(pars,model_data,scale)
#'  
#'	        scale_dm(model_data,scale)
#' 
#'          scale_par(par,scale)
#'
#'          unscale_par(par,scale)
#' 
#' @aliases set_scale scale_dm scale_par unscale_par
#' @param pars character vector of parameter names
#' @param par list of parameter vectors or vector of parameter values
#' @param scale list or vector of parameter scales
#' @param model_data list of data/design objects
#' @return List of scale values for set_scale, model.data with scaled design matrices for scale_dm,
#' vector of scaled parameter values for scale_par, and list of unscaled parameter vectors for unscale_par
#' @author Jeff Laake 
set_scale=function(pars,model_data,scale)
{
	scale.list=vector("list",length(pars))
	names(scale.list)=pars
	if(!is.null(scale)&&!is.list(scale)&&all(scale==1))
	{
		for(parx in pars)
			if(!is.null(model_data[[paste(parx,".dm",sep="")]]))
			scale.list[[parx]]=rep(1,ncol(model_data[[paste(parx,".dm",sep="")]]))
	}else
	{
		for(parx in pars)
		{
			if(is.null(scale[[parx]]))
				if(!is.null(model_data[[paste(parx,".dm",sep="")]]))
				scale.list[[parx]]=apply(model_data[[paste(parx,".dm",sep="")]],2,function(x) mean(x[x!=0]))
			else
			{
				if(length(scale[[parx]])==1)
					if(!is.null(model_data[[paste(parx,".dm",sep="")]]))
				       scale.list[[parx]]=rep(scale[[parx]],ncol(model_data[[paste(parx,".dm",sep="")]]))
			    else
				   if(length(scale[[parx]])!=ncol(model_data[[paste(parx,".dm",sep="")]]))
					   stop(paste("For",parx,"length of scale does not match length of parameters\n"))
				   else
					   scale.list[[parx]]=scale[[parx]]				   
			}		
		}
	}
	for(parx in pars)
	if(!is.null(model_data[[paste(parx,".dm",sep="")]]))
	   names(scale.list[[parx]])=colnames(model_data[[paste(parx,".dm",sep="")]])	
	return(scale.list)
}
scale_dm=function(model_data,scale)
{
	pars=names(scale)
	for(parx in pars)
		if(!is.null(model_data[[paste(parx,".dm",sep="")]]))
		model_data[[paste(parx,".dm",sep="")]]=t(t(as.matrix(model_data[[paste(parx,".dm",sep="")]]))/scale[[parx]])
    return(model_data)
}
scale_par=function(par,scale)
{
	pars=names(scale)
	for(parx in pars)
		par[[parx]]=par[[parx]]*scale[[parx]]
    return(unlist(par,use.names=FALSE))
}
unscale_par=function(par,scale)
{
	if(!is.list(par))
	{
	  pars=names(scale)
	  snames=factor(unlist(sapply(names(scale),function(x) rep(x,length(scale[[x]])))),levels=pars)
	  par.list=split(par,snames)
	} else
	{
	  pars=names(par)
	  par.list=par
	}
	for(parx in pars)
	{
		names(par.list[[parx]])=names(scale[[parx]])
		par.list[[parx]]=par.list[[parx]]/scale[[parx]]
	}
    return(par.list)
}
