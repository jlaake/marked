#' Create design data for simulation
#'  
#' Creates a dataframe that can be used with formulae to specify models for parameters in a c-r simulation experiment. 
#' 
#' A c-r data set is composed of individuals  which can be either caught or not caught on a series of occasions through time.  Each individual
#' is part of a cohort defined by the time of initial capture (or release).  Animals are assumed to age by
#' 1 unit in the interval between each capture occasion from their initial age at first capture/release.
#' If the length of initial.age =1, all animals assumed to be caught at same initial age (eg pups).  If it
#' matches length of cohorts, each cohort can have a different initial age and if it matches number of animals
#' each animal can have a different initial age. There are 1 or more records in the dataframe for each animal. 
#' For animals in the first cohort there are num.cohort records and for animals in the second cohort there are
#' num.cohorts-1 records, etc.  As an example, consider 3 cohorts with 1 animal in each cohort all starting at 
#' age 0.The dataframe would be:
#'
#' id  cohort  time  age
#'  1   1       1     0
#'  1   1       2     1
#'  1   1       3     2
#'  1   2       2     0
#'  1   2       3     1
#'  1   3       3     0
#'
#' Added covariates can be added by including another column in the dataframe.  The added covariates can be 
#'  specific to a cohort, time, age or individual.
#'  
#' @param num.cohorts number of cohorts
#' @param cohort.sizes size of each cohort; if only a single value each cohort has same size
#' @param initial.age if specified an age field is included as defined in details
#' @param num.occ number of occasions which must be as large as num.cohorts; most useful for num.cohorts=1; num.occ >1
#' @return dataframe with columns id and design data fields: cohort and time and optionally age  
#' @export create.simdesign.data
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' 
create.simdesign.data <- function(num.cohorts,cohort.sizes,initial.age=NULL,num.occ=NULL)
{
#       
# Setup cohorts and cohort sizes
#
	if(num.cohorts>1 & length(cohort.sizes)==1)
		cohort.sizes=rep(cohort.sizes,num.cohorts)
	else
	if(num.cohorts !=length(cohort.sizes))
		stop("Number of cohorts doesn't match vector of cohort sizes\n")
	if(is.null(num.occ))
		num.occ=num.cohorts
	else
	   if(num.occ<num.cohorts)stop("num.occ cannot be less than num.cohorts")
#
# Setup initial ages if not null
#
	if(!is.null(initial.age))
	{
		age.diff=FALSE
		if(num.cohorts>1 & length(initial.age)==1)
			initial.age=rep(initial.age,num.cohorts)
		else
		if(num.cohorts !=length(initial.age))
		{
			if(sum(cohort.sizes)!=length(initial.age))
				stop("Length of initial.age doesn't match number of cohorts or number of animals\n")
			else
				age.diff=TRUE
		}
	}
#
# Loop over cohorts and create design data
#
	ddl=NULL
	for (i in 1:num.cohorts)
	{
		time=as.vector(apply(as.matrix((1+i-1):num.occ),1,FUN=function(x,size) {return(rep(x,size))},size=cohort.sizes[i]))
#
# Only include age if initial.age is not null
#
		if(is.null(initial.age))
			ddl=rbind(ddl,data.frame(id=rep((sum(cohort.sizes[1:i])-cohort.sizes[i]+1):sum(cohort.sizes[1:i]),(num.occ-i+1)), cohort=rep(i, cohort.sizes[i]*(num.occ-i+1)), 
							time=time))
		else
		if(!age.diff)
			ddl=rbind(ddl,data.frame(id=rep((sum(cohort.sizes[1:i])-cohort.sizes[i]+1):sum(cohort.sizes[1:i]),(num.occ-i+1)), cohort=rep(i, cohort.sizes[i]*(num.occ-i+1)), 
							time=time,age=initial.age[i]+time-i))
		else
			ddl=rbind(ddl,data.frame(id=rep((sum(cohort.sizes[1:i])-cohort.sizes[i]+1):sum(cohort.sizes[1:i]),(num.occ-i+1)), cohort=rep(i, cohort.sizes[i]*(num.occ-i+1)), 
							time=time,age=initial.age[(sum(cohort.sizes[1:i])-cohort.sizes[i]+1):sum(cohort.sizes[1:i])]+time-i))
	}
	row.names(ddl)=1:dim(ddl)[1]
	return(ddl)
}
create.parmat <- function(x,nocc,N)
#
# create.parmat   - used in sim code to setup a matrix of parameter
#                   values with rows corresponding to simulated animals and columns
#                   are for each occasion
#
#
{
	if(is.matrix(x))
	{
		if(dim(x)[1]!=N | dim(x)[2] != (nocc-1))
			stop("Number of rows in parameter matrix is not equal to N or number of columns is not equal to nocc-1\n")
	}
	else
	{
		if(length(x)==1 | length(x)==(nocc-1) )
			x=matrix(x,nrow=N,ncol=(nocc-1),byrow=T)
		else
		if(length(x)==N)
			x=matrix(x,nrow=N,ncol=(nocc-1))
		else
			stop(paste("Invalid vector length - must be 1, nocc-1 or N\n",x))
	}
	return(x)
}


