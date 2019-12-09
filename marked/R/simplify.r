simplify_indices=function(x)
{
    #if(is(x,"dgCMatrix"))x=as.matrix(x)
    nms = colnames(x)
    df=data.table(x)
    u_rows=which(!duplicated(df))
    u_df = df[u_rows,,drop=FALSE]
    u_df$idx = 1:nrow(u_df)
    out_idx = merge(df, u_df, by=nms, sort=F)
    return(list(indices=out_idx$idx,set=u_rows))
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
			slist=simplify_indices(ddl[[parname]][,fields,drop=FALSE])
		ddl[[parname]]=ddl[[parname]][slist$set,]
		ddl[[paste(parname,"indices",sep=".")]]=slist$indices
	}
	return(ddl)
}

