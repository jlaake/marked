#' HMM computation demo functions
#' 
#' Uses dataframe and arguments to construct HMM state vectors alpha and phi for
#' demonstration.#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param id id number for calculation
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model (hmm models only) 
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional vector of initial values for beta parameters; if
#' named from previous analysis only relevant values are used
#' @param groups Vector of names factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce
#' computation
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories and design matrices; amount used is 8*chunk_size/1e6 MB (default
#' 80MB)
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata. For any unobserved strata p=0..
#' @param state.names names for states used to label output
#' @param obs.names names for observations used to lablel output 
#' @param ... optional arguments passed to js or cjs and optimx
#' @return hmm demo list
#' @author Jeff Laake
#' @export hmmDemo
#' @useDynLib marked
#' @keywords models
#' @examples
#' 
#' # cormack-jolly-seber model
#' data(dipper)
#' #note id values will not match row numbers in dipper because
#' #capture histories of "0000001" are removed. But if you process data and
#' #pass it to hmmDemo then it will match dipper.proc$data 
#' dipper.proc=process.data(dipper,model="hmmCJS")
#' x=hmmDemo(dipper.proc,id=45,state.names=c("Alive","Dead"),obs.names=c("Missed","Seen"))
#' par(mfrow=c(2,1))
#' barplot(t(x$alpha),beside=TRUE,names.arg=x$chforwardstrings)
#' barplot(t(x$phi),beside=TRUE,names.arg=x$chforwardstrings)
#' #data(mstrata)
#' #mstrata$freq=1
#' #x=hmmDemo(mstrata,id=1,model="hmmMSCJS")
#' #rowSums(x$alpha*x$beta)
#' #exp(x$lnl)
#' # state probabilities given the data
#' #x$alpha*x$beta/exp(x$lnl)
hmmDemo <- function(data,id,ddl=NULL,begin.time=1,model="hmmCJS",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
		groups = NULL, time.intervals = NULL,accumulate=FALSE,chunk_size=1e7, strata.labels=NULL,state.names=NULL,obs.names=NULL,...)
{
	model=toupper(model)
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
# 
	if(is.null(data$data))
	{
		if(!is.null(ddl))
		{
			warning("Warning: specification of ddl ignored, as data have not been processed")
			ddl=NULL
		}
		message("Model: ",model,"\n")
		message("Processing data\n")
		flush.console()
		data.proc=process.data(data,begin.time=begin.time, model=model,mixtures=1, 
				groups = groups, age.var = NULL, initial.ages = NULL, 
				time.intervals = time.intervals,nocc=NULL,accumulate=accumulate,strata.labels=strata.labels)
	}   
	else
	{
		data.proc=data
		model=data$model
	}
#
# If the design data have not been constructed, do so now
#
	if(is.null(ddl)) 
	{
		message("Creating design data.\n")
		flush.console()
		ddl=make.design.data(data.proc,design.parameters)
	} else
		design.parameters=ddl$design.parameters
#
# Setup parameter list
#
	number.of.groups=1
	if(!is.null(data.proc$group.covariates))number.of.groups=nrow(data.proc$group.covariates)
	par.list=setup.parameters(data.proc$model,check=TRUE)
#
# Check validity of parameter list; stop if not valid
#
	if(!valid.parameters(model,model.parameters)) stop()
	parameters=setup.parameters(data.proc$model,model.parameters,data$nocc,number.of.groups=number.of.groups)
	parameters=parameters[par.list]
# See if any formula contain random effects and set re
	re=FALSE
	for (i in 1:length(parameters))
	{
		if(is.null(parameters[[i]]$formula)) parameters[[i]]$formula=~1
		mlist=proc.form(parameters[[i]]$formula)
		if(!is.null(mlist$re.model))re=TRUE
	}
#  setup fixed values 
	ddl=set.fixed(ddl,parameters)
# Create design matrices for each parameter
	dml=create.dml(ddl,model.parameters=parameters,design.parameters=design.parameters,chunk_size=1e7)
# For HMM call set.initial to get ptype and set initial values
	if(substr(model,1,3)=="HMM")
		initial.list=set.initial(names(dml),dml,initial)
	else
		stop(" This function is only for Hidden Markov Models")
	if(is.null(data.proc$strata.list))		
		object=loglikelihood(unlist(initial.list$par),type=initial.list$ptype,x=data.proc$ehmat,id,m=data.proc$m,T=data.proc$nocc,start=data.proc$start,freq=data.proc$freq,
					fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,parameters=parameters)
	else
	{
		m=list(ns=length(data.proc$strata.list$states),na=length(data.proc$strata.list[[names(data.proc$strata.list)[names(data.proc$strata.list)!="states"]]]))
		object=loglikelihood(unlist(initial.list$par),type=initial.list$ptype,x=data.proc$ehmat,id,m=data.proc$m,T=data.proc$nocc,start=data.proc$start,freq=data.proc$freq,
				fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,parameters=parameters)
	}
	if(is.null(state.names))state.names=1:ncol(object$alpha)
	if(is.null(obs.names))obs.names=1:dim(object$dmat)[1]
    dimnames(object$alpha)[2]=list(state.names)
	dimnames(object$phi)=dimnames(object$alpha)
	dimnames(object$v)=dimnames(object$alpha)
	dimnames(object$gamma)[1:2]=list(state.names,state.names)
	names(dimnames(object$gamma))=c("From_state","To_state","Occasion")
	dimnames(object$dmat)[1:2]=list(obs.names,state.names)
	names(dimnames(object$dmat))=c("Observation","State","Occasion")
	colnames(object$alpha)=state.names
	colnames(object$phi)=state.names
	names(dimnames(object$alpha))=c("Occasion","State")
	names(dimnames(object$phi))=names(dimnames(object$alpha))
	names(dimnames(object$v))=names(dimnames(object$alpha))
	ch=strsplit(data.proc$data$ch[data.proc$data$id==id],",")[[1]]
	object$chfowardstrings=sapply(1:length(ch),function(x)paste(ch[1:x],collapse=""))
	object$chbackwardstrings=sapply(1:length(ch),function(x)paste(rev(ch)[1:x],collapse=""))
	return(object)
}
