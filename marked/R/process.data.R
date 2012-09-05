#' Process encounter history dataframe for MARK analysis
#' 
#' Prior to analyzing the data, this function initializes several variables
#' (e.g., number of capture occasions, time intervals) that are often specific
#' to the capture-recapture model being fitted to the data.  It also is used to
#' 1) define groups in the data that represent different levels of one or more
#' factor covariates (e.g., sex), 2) define time intervals between capture
#' occasions (if not 1), and 3) create an age structure for the data, if any.
#' 
#' For examples of \code{data}, see \code{\link{dipper}}. The structure of the
#' encounter history and the analysis depends on the analysis model to some
#' extent. Thus, it is necessary to process a dataframe with the encounter
#' history (\code{ch}) and a chosen \code{model} to define the relevant values.
#' For example, number of capture occasions (\code{nocc}) is automatically
#' computed based on the length of the encounter history (\code{ch}) in
#' \code{data}. Currently, only 2 types of models are accepted in marked: cjs and js.  
#' The default time interval is unit time (1) and if this is
#' adequate, the function will assign the appropriate length.  A processed data
#' frame can only be analyzed using the model that was specified.  The
#' \code{model} value is used by the functions \code{\link{make.design.data}}
#' and \code{\link{crm}} to define the model structure as it relates to the
#' data. Thus, if the data are going to be analysed with different underlying
#' models, create different processed data sets with the model name as an
#' extension.  For example, \code{dipper.cjs=process.data(dipper)}.
#' 
#' This function will report inconsistencies in the lengths of the capture
#' history values and when invalid entries are given in the capture history.
#' 
#' The argument \code{begin.time} specifies the time for the first capture
#' occasion.  This is used in creating the levels of the time factor variable
#' in the design data and for labelling parameters. If the \code{begin.time}
#' varies by group, enter a vector of times with one for each group. Note that
#' the time values for survivals are based on the beginning of the survival
#' interval and capture probabilities are labeled based on the time of the
#' capture occasion.  Likewise, age labels for survival are the ages at the
#' beginning times of the intervals and for capture probabilities it is the age
#' at the time of capture/recapture.
#' 
#' \code{groups} is a vector of variable names that are contained in
#' \code{data}.  Each must be a factor variable. A group is created for each
#' unique combination of the levels of the factor variables.  In the first
#' example given below \code{groups=c("sex","age","region")}. which creates
#' groups defined by the levels of \code{sex}, \code{age} and \code{region}.
#' There should be 2(sexes)*3(ages)*4(regions)=24 groups but in actuality there
#' are only 16 in the data because there are only 2 age groups for each sex.
#' Age group 1 and 2 for M and age groups 2 and 3 for F.  This was done to
#' demonstrate that the code will only use groups that have 1 or more capture
#' histories unless \code{allgroups=TRUE}.
#' 
#' The argument \code{age.var=2} specifies that the second grouping variable in
#' \code{groups} represents an age variable.  It could have been named
#' something different than age. If a variable in \code{groups} is named
#' \code{age} then it is not necessary to specify \code{age.var}.
#' \code{initial.age} specifies that the age at first capture of the age levels
#' is 0,1 and 2 while the age classes were designated as 1,2,3. The actual ages
#' for the age classes do not have to be sequential or ordered, but ordering
#' will cause less confusion.  Thus levels 1,2,3 could represent initial ages
#' of 0,4,6 or 6,0,4. The default for \code{initial.age}
#' is 0 for each group, in which case, \code{age} represents time since marking
#' (first capture) rather than the actual age of the animal.
#' 
#' The following variable names are reserved and should be used as follows:
#' id (animal id)
#' ch(capture history)
#' freq (number of animals with that ch/data)
#' The following variable names are reserved and should not be used in the data:
#' age,time,cohort,Age,Time,Cohort,Y,Z
#' 
#' @aliases process.data accumulate_data
#' @usage 	process.data(data,begin.time=1,model="cjs",mixtures=1,groups=NULL,allgroups=FALSE,age.var=NULL,
#'               initial.ages=c(0),time.intervals=NULL,nocc=NULL,accumulate=TRUE)
#' 
#' @param data A data frame with at least one field named \code{ch} which is
#' the capture (encounter) history stored as a character string. \code{data}
#' can also have a field \code{freq} which is the number of animals with that
#' capture history. The default structure is freq=1 and it need not be included
#' in the dataframe. \code{data} can also contain an arbitrary number of
#' covariates specific to animals with that capture history.
#' @param begin.time Time of first capture occasion or vector of times if
#' different for each group
#' @param model Type of analysis model. See \code{mark in RMark} for a list of
#' possible values for \code{model}
#' @param mixtures Number of mixtures in closed capture models with
#' heterogeneity
#' @param groups Vector of factor variable names (in double quotes) in
#' \code{data} that will be used to create groups in the data. A group is
#' created for each unique combination of the levels of the factor variables in
#' the list.
#' @param allgroups Logical variable; if TRUE, all groups are created from
#' factors defined in \code{groups} even if there are no observations in the
#' group
#' @param age.var An index in vector \code{groups} for a variable (if any) for
#' age
#' @param initial.ages A vector of initial ages that contains a value for each
#' level of the age variable \code{groups[age.var]}
#' @param time.intervals Vector of lengths of time between capture occasions
#' @param nocc number of occasions for Nest type; either nocc or time.intervals
#' must be specified
#' @param accumulate if TRUE, aggregates data with same values and creates freq field for count of records
#' @return from \code{process.data} processed.data (a list with the following elements)
#' \item{data}{original raw dataframe with group factor variable added if
#' groups were defined} \item{model}{type of analysis model (eg, "cjs" or "js")}
#' \item{freq}{a dataframe of frequencies (same # of rows
#' as data, number of columns is the number of groups in the data. The column
#' names are the group labels representing the unique groups that have one or
#' more capture histories.} \item{nocc}{number of capture occasions}
#' \item{time.intervals}{length of time intervals between capture occasions}
#' \item{begin.time}{time of first capture occasion} \item{initial.ages}{an initial age for
#' each group in the data; Note that this is not the original argument but is a
#' vector with the initial age for each group. In the first example below
#' \code{proc.example.data$initial.ages} is a vector with 16 elements as
#' follows 0 1 1 2 0 1 1 2 0 1 1 2 0 1 1 2} \item{group.covariates}{factor covariates used to define groups}
#' from accumulate_data a dataframe with same column structure as argument with addition of freq (if not any)
#' and reduced to unique rows with freq accumulating number of records. 
#' @author Jeff Laake
#' @export process.data accumulate_data
#' @seealso \code{\link{dipper}},\code{\link{crm}}
#' @keywords utility
#' @examples
#' 
#' 
#' data(dipper)
#' dipper.process=process.data(dipper)
#' accumulate_data(dipper)
#' 
accumulate_data <- function(data)
{
	x <- data[,names(data)!="freq"]
	nx <- nrow(x)
	if(is.null(data$freq))data$freq=rep(1,nrow(data))
	pasted.data=apply(x, 1, paste, collapse = "")
	freq=sapply(split(data$freq, pasted.data),sum)
	x=unique(x[order(pasted.data),])
	x$freq=freq
	cat(nx,"capture histories collapsed into ",nrow(x),"\n")
	return(x)	
}
process.data <-
function(data,begin.time=1,model="CJS",mixtures=1,groups=NULL,allgroups=FALSE,age.var=NULL,
initial.ages=c(0),time.intervals=NULL,nocc=NULL,accumulate=TRUE)
{
   if(model%in%c("cjs","js"))model=toupper(model)
   dataname=substitute(data)
#
#  Compute number of occasions and check validity of model
#
   if(is.null(data$ch))
     stop("Field ch is missing in ",substitute(data))
   ch.lengths=nchar(data$ch)
   nocc=median(ch.lengths)
   if(any(ch.lengths!=nocc))
   {
        stop(paste("\nCapture history length is not constant. ch must be a character string",
            "\n row numbers with incorrect ch length",paste(row.names(data[ch.lengths!=nocc,]),collapse=","),"\n"))
   }
#
#  Setup model
#
   model.list=setup.model(model,nocc,mixtures)
   ch.values=unique(unlist(strsplit(data$ch,"")))
   if(any(!ch.values%in%c("0","1",".")))
      stop(paste("\nIncorrect ch values in data:",paste(ch.values,collapse=""),"\n",sep=""))
   nocc=model.list$nocc
   nocc.secondary=NULL
   num=model.list$num
#
#     If time intervals specified make sure there are nocc-1 of them
#     If none specified assume they are 1
#
   if(is.null(time.intervals))
      time.intervals=rep(1,nocc+model.list$num)
   else
      if(length(time.intervals)!=(nocc+num))
          stop("Incorrect number of time intervals")
   mixtures=model.list$mixtures
#
#  Get number of factors to create groups
#
   if(is.null(groups))
     number.of.factors=0
   else
     number.of.factors=length(groups)
#
# Accumulate non-unique data; if accumulate=F and dataframe has freq>1 expand and add id field
# if null
# 
if(!is.null(data$Freq)) names(data)[which("Freq"== names(data))]="freq"
if(model=="probitCJS") accumulate=F
if(accumulate)
	data=accumulate_data(data)
else
{
	if(is.null(data$freq))
		data$freq=1
	else
	{
		data=data[rep(1:nrow(data),times=data$freq),]
		data$freq=1
	}
	if(is.null(data$id))
		data$id=1:nrow(data)
	else
		data=data[order(data$id),]
}
#
#  Get number of records in data set
#
number.of.ch=nrow(data)
#
#  If there are no factors then
#     if already has freq variable return the input data set as a list
#     otherwise add the freq variable with each value = 1 and return as a list
#  If model=js, then add dummy data for non-captured 
#
if(number.of.factors==0)
{
    if(model=="js")
    {
       data=add.dummy.data(data,nocc=nocc,group.covariates=NULL)     
       number.of.ch=nrow(data)
    }
    return(list(data=data,model=model,mixtures=mixtures,
                   freq=matrix(data$freq,ncol=1,dimnames=list(1:number.of.ch,"group1")),
                   nocc=nocc, nocc.secondary=nocc.secondary,time.intervals=time.intervals,begin.time=begin.time,
                   initial.ages=initial.ages[1],group.covariates=NULL))
}
#
#   If there are one or more in the group factor list then
#     make sure each is a factor variable in the data set and compute number
#         of levels for each factor, cumlevels and factor matrix
#     if not a factor variable - stop with error message
# 
else
{
  number.of.groups=1
  n.levels=rep(0,number.of.factors)
  facmat=NULL
  faclabs=list()
  for (i in 1:number.of.factors)
  {
    vari=data[,groups[i]]
    if(!is.factor(vari))
        stop(paste("\n ",groups[i]," is not a factor variable\n"))
     else
     {
        n.levels[i]=length(levels(vari))
        facmat=cbind(facmat,as.numeric(vari)-1)
        faclabs[[i]]=levels(vari)
     }       
  }
  cumlevels=cumprod(n.levels)
  number.of.groups=cumlevels[length(cumlevels)]

#  If age.var is specified, make sure it is valid and that the number of 
#  initial.ages matches number of levels of identified variable
#
   if(is.null(age.var))
      age.var=match("age",groups)
   if(!is.na(age.var))
   {
      if(age.var>length(groups) | age.var<1)
         stop("Invalid age variable. Must be between 1 and ",length(groups))
      if(is.null(initial.ages))
         stop("initial.ages must be specified if age.var is specified")
      else
         if(!is.numeric(initial.ages) | (length(initial.ages)!=n.levels[age.var] &length(initial.ages)>1) )
           stop(paste("intial.ages must be numeric and match length of levels of",groups[age.var]))
   }
#
#  Next compute the group number for each capture history
#
   if(number.of.factors==1)
      data$group=facmat+1
   else
      if(number.of.factors==2)
         data$group=facmat[,2]*cumlevels[1]+facmat[,1]+1
      else
         data$group=facmat[,2:number.of.factors]%*%cumlevels[1:(number.of.factors-1)]+facmat[,1]+1
#
#  Next create frequency matrix for groups   
#
  freqmat=matrix(0,nrow=number.of.ch,ncol=number.of.groups)
  for(i in 1:number.of.ch)
  {
     freqmat[i,data$group[i]]=data$freq[i]
  }
#
#  If allgroups=FALSE, recompute number of groups and group number based on groups with 1 or more capture histories
#
  if(!allgroups)
  {
     test.freq=freqmat
     test.freq[test.freq!=0]=1
     counts = apply(test.freq, 2, sum)
     newgroups=rep(0,number.of.groups)
     index=1
     for (i in 1:number.of.groups)
        if(counts[i]>0)
        {
           newgroups[i]=index
           index=index+1
        }     
     data$group=as.factor(newgroups[data$group])
     freqmat=freqmat[,counts>0]
     number.of.groups=index-1
  }
#
#  Check to make sure length of begin.time is either 1 or equal to the
#  number of groups
#
  if(length(begin.time)!=1 & length(begin.time)!=number.of.groups)
    stop("length of begin.time must either be 1 or match number of groups")
#
#  Create group labels
#  
  labs=expand.grid(faclabs)
  if(!allgroups)labs=as.matrix(labs[counts>0,])
#
#  If age.var has not been set, initial ages are set to 0
#
  if(is.na(age.var))
    init.ages=rep(initial.ages[1],number.of.groups)
  else
  {
    if(length(initial.ages)==1)
       initial.ages=rep(initial.ages,length(levels(as.factor(labs[,age.var]))))
    init.ages = initial.ages[as.numeric(factor(labs[,age.var],levels=unique(faclabs[[age.var]])))]
  }
  grouplabs=rep(" ",number.of.groups)
  for (i in 1:number.of.groups)
     grouplabs[i]=paste(groups,labs[i,],sep="",collapse=".") 
  freqmat=as.data.frame(freqmat)
  names(freqmat)=grouplabs
#
#  Store labs as group covariates; set levels to the same as in data
#  
  group.covariates=as.data.frame(labs)
  names(group.covariates)=groups
  for (i in 1:dim(group.covariates)[2])
     group.covariates[,i]=factor(group.covariates[,i],levels=levels(data[,groups[i]]))
#
# Return data as a list with original dataframe and frequency matrix
#
  if(model=="js")
     data=add.dummy.data(data,nocc,group.covariates)     
  data$initial.age=init.ages[data$group]
  return(list(data=data,model=model,mixtures=mixtures,freq=freqmat,
                   nocc=nocc, nocc.secondary=nocc.secondary, time.intervals=time.intervals,begin.time=begin.time,
                   initial.ages=init.ages,group.covariates=group.covariates))
}
}
add.dummy.data=function(data,nocc,group.covariates)
{
	if(is.null(group.covariates))
		number.of.groups=1
	else
		number.of.groups=nrow(group.covariates)
	if(!is.null(group.covariates))
	{
		xlist=split(data,data[,names(group.covariates),drop=FALSE])
		xlist=xlist[as.vector(sapply(xlist,function(x) nrow(x)))>0]
	} else
	{                              
		xlist=data
	}
	numvar=sapply(data,is.numeric)
	numvar=numvar[names(data)!="freq"]
	if(any(numvar))
	{
		numvar=names(data[,names(data)!="freq"])[numvar]
		if(!is.null(group.covariates))
		{
			xmeans=sapply(xlist,function(x) sapply(subset(x,select=numvar),mean))
			if(length(numvar)==1)
			{
				dd=data.frame(group=1:nrow(group.covariates))
				dd[,numvar]=xmeans
			}else
			{
				dd=data.frame(group=1:nrow(group.covariates))
				dd=cbind(dd,t(xmeans))
			}
			xx=merge(cbind(data.frame(group=1:nrow(group.covariates),group.covariates)),
					dd,by.x="group",all.x=TRUE)
		}
		else
		{
			xmeans=colMeans(subset(data,select=numvar))    
			if(ncol(t(xmeans))==0)
				xx=data.frame(N=1)
			else
				xx=data.frame(xmeans)  
		}
	}else
	{
		numvar=NULL
		if(!is.null(group.covariates))
			xx=cbind(group.covariates,group=factor(1:nrow(group.covariates)))
		else
			xx=NULL
	}
	chmat=matrix(0,nrow=nocc,ncol=nocc)
	diag(chmat)=1
	ch=apply(chmat,1,paste,collapse="")
	ch=rep(ch,number.of.groups)
	if(is.null(data$freq))
		data$freq=1
	if(!is.null(xx))
	{
		if(number.of.groups==1)
		{
			data=subset(data,select=c("ch","freq",numvar,names(group.covariates)))
			dummy.data=cbind(data.frame(ch=ch,freq=rep(0,length(ch))),matrix(t(xx),byrow=T,ncol=length(numvar),nrow=length(ch)))
			names(dummy.data)=c("ch","freq",numvar)             
		}else
		{
			data=subset(data,select=c("ch","freq","group",numvar[numvar!="group"],names(group.covariates)))
			xx=subset(xx,select=!names(xx)%in%"group")
			dummy.data=cbind(data.frame(ch=ch,freq=rep(0,length(ch))),group=factor(rep(1:number.of.groups,each=nocc)),
					xx[rep(1:number.of.groups,each=nocc),,drop=FALSE])
			names(dummy.data)=c("ch","freq","group",names(group.covariates),numvar[numvar!="group"])             
		}
	}
	else
	{
		data=subset(data,select=c("ch","freq"))
		dummy.data=data.frame(ch=ch,freq=rep(0,length(ch)))
		names(dummy.data)=c("ch","freq")             
	}
	row.names(dummy.data)=NULL
	return(rbind(data,dummy.data))
}
