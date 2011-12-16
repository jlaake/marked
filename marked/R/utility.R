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
#' @param x processed data list - result from process.data in marked
#' @param object model fitted by marked
#' @param age at which Phi or p should be shown across time
#' @param time at which Phi or p should be shown across ages
#' @export resight.matrix 
#' @export naive.survival
#' @export Phi.mean
#' @export p.mean
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
Phi.mean=function(object,age=0,time=NULL)
{	
	if(is.null(object$real$sex))object$real$sex="All"
	if(is.null(time))
	{	
		with(object$real[object$real$Time>=object$real$Cohort&object$real$Age==age,],tapply(Phi,list(sex,Phi.time),mean))
	} else
	{
		with(object$real[object$real$Phi.time==time,],tapply(Phi,list(sex,Phi.age),mean))
	}
}
p.mean=function(object,age=0,time=NULL)
{	
	if(is.null(object$real$sex))object$real$sex="All"
	if(is.null(time))
	{
		with(object$real[object$real$Time>=object$real$Cohort&object$real$Age==age,],tapply(p,list(sex,p.time),mean))
	} else
	{	
		with(object$real[object$real$p.time==time,],tapply(p,list(sex,p.age),mean))
	}
}
Phi.boxplot=function(object,age=0,time=NULL,sex=NULL){
	if(!is.null(sex))
	{
	   if(is.null(time))
	   {
		   boxplot(Phi~Phi.time,data=object$real[object$real$Age==age&object$real$sex==sex,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age,"and sex ",sex))
	   } else
	   {
		   boxplot(Phi~Phi.age,data=object$real[object$real$Phi.time==time&object$real$sex==sex,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time,"and sex ",sex))
		   
	   }
   }else
   {
	   if(is.null(time))
	   {
		   boxplot(Phi~Phi.time,data=object$real[object$real$Age==age,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age))
	   } else
	   {
		   boxplot(Phi~Phi.age,data=object$real[object$real$Phi.time==time,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time))
	   }
	}
}

p.boxplot=function(object,age=0,time=NULL,sex=NULL){
	if(!is.null(sex))
	{
		if(is.null(time))
		{
			boxplot(p~p.time,data=object$real[object$real$Age==age&object$real$sex==sex,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age,"and sex ",sex))
		} else
		{
			boxplot(p~p.age,data=object$real[object$real$p.time==time&object$real$sex==sex,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time,"and sex ",sex))
			
		}
	}else
	{
		if(is.null(time))
		{
			boxplot(p~p.time,data=object$real[object$real$Age==age,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age))
		} else
		{
			boxplot(p~p.age,data=object$real[object$real$p.time==time,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time))
		}
	}
}

