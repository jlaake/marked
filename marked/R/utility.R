#' Various utility functions
#' 
#' Several functions have been added to help visualize data including resight.matrix
#' which provides for each cohort the number of releases on the diagonal and the number 
#' resighted on each occasion in the upper-triangular matrix. naive.survival
#' provides a naive survival estimate at each time for each cohort.  The estimate for
#' time i is the number resighted at time i+1 or later divided by the number seen at time i or later.
#' These values should be interpreted cautiously because they are influenced by capture probability in
#' year i but it is useful to identify particularly high or low survival values. Functions Phi.mean
#' and p.mean compute average real parameter values Phi or p across time for a single age or across ages for a
#' single time.
#' 
#' @usage resight.matrix(x)
#'           naive.survival(x,...)
#' @param x processed data list - result from process.data in marked or real estimates from fitted model
#' @param age at which Phi or p should be shown across time
#' @param time at which Phi or p should be shown across ages
#' @export resight.matrix 
#' @export naive.survival
#' @export Phi.mean
#' @export p.mean
#' @export function.wrapper
#' @export fx.aic
#' @export fx.par.count
#' @return matrix of values with cohort and year labels
#' @author Jeff Laake
#' @keywords utility
resight.matrix=function(x)
{
	zz=strsplit(x$data$ch,"")
	zz=t(sapply(zz,as.numeric))
	resight.matrix=t(sapply(split(as.data.frame(zz),x$data$cohort),colSums))
	colnames(resight.matrix)=c(0,cumsum(x$time.intervals))+x$begin.time
	return(resight.matrix)
}
naive.survival=function(x)
{
	ll=process.ch(x$data$ch,all=TRUE)
	surv.matrix=t(sapply(split(as.data.frame(ll$First-ll$Lplus),x$data$cohort),colSums))
	naive.S=matrix(NA,nrow=nrow(surv.matrix),ncol=ncol(surv.matrix))
	naive.S=surv.matrix[,2:(ncol(surv.matrix))]/surv.matrix[,1:(ncol(surv.matrix)-1)]
	naive.S[is.infinite(naive.S)]=NA
	naive.S[is.nan(naive.S)]=NA
	rownames(naive.S)=levels(x$data$cohort)
	colnames(naive.S)=x$begin.time+ c(0,cumsum(x$time.intervals[-length(x$time.intervals)]))
	return(naive.S)	
}
Phi.mean=function(x,age=0,time=NULL,age.bins=NULL,age.levels=NULL)
{	
	if(is.null(x$sex))x$sex="All"
	if(is.null(time))
	{	
		with(x[x$Time>=x$Cohort&x$Age%in%age,],tapply(Phi,list(sex,Phi.time),mean))
	} else
	{
		x$age=cut(as.numeric(x$Phi.age),age.bins,right=FALSE)
		levels(x$age)=age.levels
		with(x[x$Phi.time%in%time,],tapply(Phi,list(sex,age),mean))
	}
}
p.mean=function(x,age=0,time=NULL,age.bins=NULL,age.levels=NULL)
{	
	if(is.null(x$sex))x$sex="All"
	if(is.null(time))
	{
		with(x[x$Time>=x$Cohort&x$Age%in%age,],tapply(p,list(sex,p.time),mean))
	} else
	{	
		if(is.null(age.bins))
		   with(x[x$p.time%in%time,],tapply(p,list(sex,p.age),mean))
	    else
		{
			x$age=cut(as.numeric(x$p.age),age.bins,right=FALSE)
			levels(x$age)=age.levels
			with(x[x$p.time%in%time,],tapply(p,list(sex,age),mean))
		}
	}
}


Phi.boxplot=function(x,age=0,time=NULL,sex=NULL){
	if(!is.null(sex))
	{
	   if(is.null(time))
	   {
		   boxplot(Phi~Phi.time,data=x[x$Age==age&x$sex==sex,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age,"and sex ",sex))
	   } else
	   {
		   boxplot(Phi~Phi.age,data=x[x$Phi.time%in%time&x$sex==sex,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time,"and sex ",sex))
		   
	   }
   }else
   {
	   if(is.null(time))
	   {
		   boxplot(Phi~Phi.time,data=x[x$Age==age,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age))
	   } else
	   {
		   boxplot(Phi~Phi.age,data=x[x$Phi.time==time,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time))
	   }
	}
}

p.boxplot=function(x,age=0,time=NULL,sex=NULL){
	if(!is.null(sex))
	{
		if(is.null(time))
		{
			boxplot(p~p.time,data=x[x$Age==age&x$sex==sex,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age,"and sex ",sex))
		} else
		{
			boxplot(p~p.age,data=x[x$p.time==time&x$sex==sex,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time,"and sex ",sex))
			
		}
	}else
	{
		if(is.null(time))
		{
			boxplot(p~p.time,data=x[x$Age==age,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age))
		} else
		{
			boxplot(p~p.age,data=x[x$p.time==time,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time))
		}
	}
}

function.wrapper <- function(x,fx,base="",...)
{
	model.name=paste(x,collapse=".")
	e1=new.env()
	dots=list(...)
	sapply(1:length(dots),function(x) assign(names(dots)[x],dots[[x]],envir=e1))
	environment(fx)=e1
	eval(parse(text=paste('load(file="',base,model.name,'.rda")',sep="")),envir=e1)
	eval(parse(text=paste("result<-fx(",model.name,")",sep="")),envir=e1)	
	return(get("result",envir=e1))
}

fx.aic=function(x) x$neg2lnl/chat + 2*length(x$beta)
fx.par.count=function(x) length(grep(paste("^",par,":",sep=""),names(x$beta)))



