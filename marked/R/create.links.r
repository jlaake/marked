create.links=function(dm)
{
	links=rep(0,nrow(dm))
	possible.sin=apply(dm,1,function(x){
				w=which(x==1)
				ifelse(length(w)==1,w,NA)
			})								 
	valid.columns=(1:ncol(dm))[sapply(1:ncol(dm),function(x) all(dm[dm[,x]==1,-x]==0))]
	if(length(valid.columns)>0)
		for (i in valid.columns)
			links[possible.sin==i]=1
	return(links)
}
