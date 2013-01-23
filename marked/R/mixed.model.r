#' Mixed effect model contstruction
#' 
#' Functions that develop structures needed for a mixed effect model
#' 
#' mixed.model.admb - creates design matrices and supporting index matrices
#' for use of mixed model in ADMB
#' 
#' mixed.model - creates design matrices and supporting index matrices
#' in an alternate list format that is not as easily used in ADMB
#' 
#' mixed.model.dat - writes to data file (con) for fixed and random effect stuctures
#' 
#' @usage mixed.model.admb(f,data)
#'           mixed.model(f,data)
#'           mixed.model.dat(x,con)
#' @param f formula for mixed effect mode in the form used in lme4; ~fixed +(re1|g1) +...+(ren|gn)
#' @param data dataframe used to construct the design matrices from the formula
#' @param x list structure created by mixed.model.admb
#' @param con connection to data file which contents will be appended
#' @return mixed.model.admb returns a list with elements fixed.dm, the design matrix for
#' the fixed effects; re.dm, a combined design matrix for all of the random effects; and 
#' re.indices, matrix of indices into a single vector of random effects to be applied to the 
#' design matrix location.
#' mixed.model returns similar quantities for the random effects in a list structure (re.list) except that the indices are limited 
#' to the particular random effect grouping. May be more useful with R than ADMB.
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' 
mixed.model.admb=function(formula,data)
{
# parse formula for fixed and random effects
  mlist=proc.form(formula)
# construct design matrix for fixed effects
  fixed.dm=model.matrix(as.formula(mlist$fix.model),data)
# remainder of code is for random effects unless NULL
  reindex=0
  re.dm=NULL
  re.indices=NULL
  if(!is.null(mlist$re.model))
  {
#    Loop over each random effect component
     for(i in 1:length(mlist$re.model))
	 {
#      Make sure each variable used to define random effect group is a factor variable
  	   if(!all(sapply(model.frame(as.formula(mlist$re.model[[i]]$sub),data),is.factor)))
	       stop(paste("\n one or more variables in",mlist$re.model[[i]]$sub,"is not a factor variable\n"))
       else
	   {
#         Compute design matrix for grouping variables
          zz=model.matrix(as.formula(mlist$re.model[[i]]$sub),data)
#         Not all combinations of factor variable(s) may be used so only use those observed in the data
		  used.columns=which(colSums(zz)>0)
		  nre=length(used.columns)
#         Compute the indices for this particular grouping structure and reindex if any missing
          indices=rowSums(zz*col(zz))
		  if(nre!=ncol(zz))indices=match(indices,used.columns)
#         Compute the design matrix for the random effect formula
          zz=model.matrix(as.formula(mlist$re.model[[i]]$model),data)
#         Now shift indices to refer to a single vector of random effects across all re groupings
		  ng=max(indices)
          indices=matrix(indices,nrow=length(indices),ncol=ncol(zz))+reindex
          indices=t(t(indices)+cumsum(c(0,rep(ng,ncol(zz)-1))))
          reindex=max(indices)	
#         Bind random effect design matrices (re.dm), indices into random effects vector (re.indices) and
#         index for the random effect sigma parameter (re.sigma)
		  re.dm=cbind(re.dm,zz)
		  re.indices=cbind(re.indices,indices)
	   }
	 }
   }
   return(list(fixed.dm=fixed.dm,re.dm=re.dm,re.indices=re.indices))
}
mixed.model=function(formula,data)
{
  mlist=proc.form(formula)
  fixed.dm=model.matrix(as.formula(mlist$fix.model),data)
  re.list=NULL
  if(!is.null(mlist$re.model))
  {
     re.list=vector("list",length=length(mlist$re.model))
     for(i in 1:length(mlist$re.model))
 {
     if(!all(sapply(model.frame(as.formula(mlist$re.model[[i]]$sub),data),is.factor)))
       stop(paste("\n one or more variables in",mlist$re.model[[i]]$sub,"is not a factor variable\n"))
       else
   {
          zz=model.matrix(as.formula(mlist$re.model[[i]]$sub),data)
  used.columns=which(colSums(zz)>0)
  nre=length(used.columns)
          indices=rowSums(zz*col(zz))
  if(nre!=ncol(zz))indices=match(indices,used.columns)
          re.dm=model.matrix(as.formula(mlist$re.model[[i]]$model),data)   
      re.list[[i]]=list(indices=indices,re.dm=re.dm)
   }
 }
   }
   return(list(fixed.dm=fixed.dm,re.list=re.list))
}
mixed.model.dat=function(x,con)
{
	# number of columns of fixed dm
	write(ncol(x$fixed.dm),con,append=TRUE)
	# fixed dm
	write(t(x$fixed.dm),con,ncolumns=ncol(x$fixed.dm),append=TRUE)
	if(!is.null(x$re.dm))
	{
		# number of random effects
		write(max(x$re.indices),con,append=TRUE)
		# number of columns of re dm
		write(ncol(x$re.dm),con,append=TRUE)
		# re dm
		write(t(x$re.dm),con,ncolumns=ncol(x$re.dm),append=TRUE)
		# re indices
		write(t(x$re.indices),con,ncolumns=ncol(x$re.indices),append=TRUE)
	}
	else
	{
		# 0 no re
		write(0,con,append=TRUE)
		# number of re =0
		write(0,con,append=TRUE)
		# number of columns of re dm=0
		write(0,con,append=TRUE)
	}
	invisible()	
}



