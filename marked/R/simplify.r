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
setup_re=function(ddl,formula)
{
mixed=mixed.model.admb(formula,ddl)
nsigma=0
if(!is.null(mixed$re.dm))nsigma=ncol(mixed$re.dm)
if(!is.null(mixed$re.dm))
{
  mixed$re.indices[ddl$Time<ddl$Cohort,]=NA
  mixed=reindex(mixed,ddl$id)
  # random effect data
  krand=ncol(mixed$re.dm)
  randDM=mixed$re.dm
  randIndex=mixed$re.indices
  counts=mixed$index.counts
  mx=max(mixed$index.counts)
  idIndex=t(sapply(mixed$used.indices,function(x) return(c(x,rep(0,mx-length(x))))))
  nre=max(idIndex)
  if(nrow(idIndex)==1)idIndex=t(idIndex)
  if(krand==1)randIndex=matrix(as.vector(t(randIndex)),ncol=1)
  #   simplify random indices and dm
  s1=simplify_indices(idIndex)
  idIndex=idIndex[s1$set,,drop=FALSE]
  idIndex_i=s1$indices
  s1=simplify_indices(randIndex)
  randIndex=randIndex[s1$set,,drop=FALSE]
  randIndex_i=s1$indices
  s1=simplify_indices(randDM)
  randDM=randDM[s1$set,,drop=FALSE]
  randDM_i=s1$indices
} else {
  nre=0
  krand=0
  randDM=matrix(0,nrow=0,ncol=0)
  randIndex=matrix(0,nrow=0,ncol=0)
  counts=vector("integer",length=0)
  idIndex=matrix(0,nrow=0,ncol=0)
  idIndex_i=vector("integer",length=0)
  randDM_i=vector("integer",length=0)
  randIndex_i=vector("integer",length=0)
}
return(list(nsigma=nsigma,nre=nre,krand=krand,counts=counts,randDM=randDM,randDM_i=randDM_i,randIndex=randIndex,randIndex_i=randIndex_i,idIndex=idIndex,idIndex_i=idIndex_i))
}
