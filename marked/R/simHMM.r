#' Simulates data from Hidden Markov Model 
#' 
#' Creates a set of data from a specified HMM for capture-recapture data.
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model ("cjs" or "js" at present)
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional list (by parameter type) of initial values for beta parameters (e.g., initial=list(Phi=0.3,p=-2)
#' @param groups Vector of names of factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce computation
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata. For any unobserved strata p=0..
#' @export
#' @return dataframe with simulated data
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @examples 
#' # simulate phi(.) p(.) with 1000 Females and 100 males, 3 occasions all released on first occasion
#' df=simHMM(data.frame(ch=rep("100",2),sex=factor(c("F","M")),freq=c(1000,100),stringsAsFactors=F))
#' 
simHMM=function(data,ddl=NULL,begin.time=1,model="hmmCJS",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
		groups = NULL, time.intervals = NULL,accumulate=TRUE,strata.labels=NULL)
{ 
	setup=fitHMM(data=data,ddl=ddl,begin.time=begin.time,model=model,title=title,model.parameters=model.parameters,
			design.parameters=design.parameters,initial=initial,groups = groups, time.intervals = time.intervals, accumulate=accumulate,
			run=FALSE,strata.labels=strata.labels)
	parlist=split(setup$par,setup$ptype)
	T=setup$data$nocc
	m=setup$data$m
	ch=NULL
	df2=NULL
	for (id in as.numeric(setup$data$data$id))
	{
		dmat=setup$data$fct_dmat(id,ddl=setup$ddl,parlist,parameters=setup$model.parameters)
		gamma=setup$data$fct_gamma(id,ddl=setup$ddl,parlist,parameters=setup$model.parameters)
		# set up state with freq rows
		history=matrix(0,nrow=setup$data$data$freq[id],ncol=T)
		state=matrix(0,nrow=setup$data$data$freq[id],ncol=T)
		state[,setup$start[id,2]]=setup$start[id,setup$start[id,1]]
		history[,setup$start[id,2]]=setup$start[id,setup$start[id,1]]
		for(j in setup$start[id,2]:(T-1))
		{
			for(k in 1:m)
			{
				instate=sum(state[,j]==k)
				if(instate>0)
					state[state[,j]==k,j+1]= apply(rmultinom(instate,1,gamma[[j]][k,]),2,function(x) which(x==1))
			} 
			# use dmat to create observed sequence
			for(k in 1:m)
			{
				instate=sum(state[,j+1]==k)
				if(instate>0)
					history[state[,j+1]==k,j+1]= setup$data$ObsLevels[apply(rmultinom(instate,1,dmat[[j]][,k]),2,function(x) which(x==1))]
			}
		}
	    ch=c(ch,apply(history,1,paste,collapse=""))	
		if(is.null(df2))
			df2=setup$data$data[rep(id,setup$data$data$freq[id]),-which(names(setup$data$data)%in%c("ch","freq","id")),drop=FALSE]
		else
			df2=rbind(df2,setup$data$data[rep(id,setup$data$data$freq[id]),-which(names(setup$data$data)%in%c("ch","freq","id")),drop=FALSE])		
	}   
	df=data.frame(ch=ch)
	if(nrow(df2)==0)
	   return(df)
    else
	   return(cbind(df,df2))
}
