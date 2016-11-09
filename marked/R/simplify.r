"realign.pims" <-function(dm){
#  Value:
#
#  new.indices - a vector of new indices for the old PIM values.  The old
#                PIM values are 1:length(new.indices) and the new index is
#                the corresponding value.  For example, new.indices=c(1,1,2,2)
#                renumbers the structure 1,2,3,4 such that 1,2 are now 1
#                and 3,4 are now 2.
#
#  Get all the rows in the design matrix and paste all the values
#  in each row together.
#
	allvals=apply(dm,1,paste,collapse="")
#
#  Get all the unique rows in the design matrix and paste all the values
#  in each row together.
#
	uniquevals=unique(allvals)
#
#  Find the corresponding sets of indices by matching allvals into uniquevals
#
	new.indices=match(allvals, uniquevals)
	

	return(new.indices)
}
simplify_indices=function(x)
{
	uniquevals=apply(unique(x),1,paste,collapse="")
	allvals=apply(x,1,paste,collapse="")
	new.indices=match(allvals, uniquevals)
	return(list(indices=new.indices,set=which(!duplicated(allvals))))
}
simplify_ddl=function(ddl,parameters)
{
	for (parname in names(parameters))
	{
		fields=all.vars(parameters[[parname]]$formula)
		for(i in seq_along(fields))
			if(!fields[i]%in%colnames(ddl[[parname]])) stop(paste("\n",fields[i]," variable used in formula for",parname,", not found in data\n"))
		if(!is.null(ddl[[parname]]$fix)) fields=c(fields,"fix")
		if(length(fields)==0)
			slist=list(indices=rep(1,nrow(ddl[[parname]])),set=1)
		else
			slist=simplify_indices(cbind(ddl[[parname]][,fields]))
		ddl[[paste(parname,"id",sep=".")]]=ddl[[parname]]$id
		ddl[[parname]]=ddl[[parname]][slist$set,]
		ddl[[paste(parname,"indices",sep=".")]]=slist$indices
	}
	return(ddl)
}

