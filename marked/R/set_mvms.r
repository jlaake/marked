#' Multivariate Multistate (mvms) Specification
#' 
#' Creates list data structure from mvms specification 
#' 
#' Accepts a mvms specification which is a list with named character vectors.
#' The length of the list is the multivariate dimension. The name of each 
#' dimension is the list element name.  Each character vector specifies the 
#' one character labels for the states of that dimension and optionally a 
#' reserved character "u" to specify that there is state uncertainty for that
#' dimension. 
#' 
#' The code tests to make sure that the input mvms specification is of the correct 
#' structure and it stops with an error message if not. The code returns a list
#' structure with a number of elements described under return value below.
#' 
#' @param x a multivariate multistate (mvms) specification as described above
#' @return a list with the following elements: 1) mvms - the input specification,
#' 2) nd - the number of dimensions, 3) df - the dataframe containing all combinations 
#' of observations across dimensions including uncertain states, 4) df.states - the dataframe with all
#' combinations of states across dimensions, 5) uncertain - boolean vector with nd elements
#' indicating whether there is uncertainy in states for each dimension.
#' @author Jeff Laake
#' @export
#' @examples
#' set_mvms(list(location=c("A","B","C"),repro_status=c("N","P","u")))
set_mvms=function(x)
{
	if(!is.list(x))
		stop("mvms specification must be a list")
	if(any(names(x)==""))
	    stop("all dimensions of mvms specification must be named")
	if(any(!sapply(x,function(x) is.vector(x))))
		stop("each set of states must be specified as a vector")
	if(any(!sapply(x,function(x) mode(x)=="character")))
	    stop("each set of states must be specified as a character vector")
	nd=length(x)
	uncertain=sapply(x,function(x)"u"%in%x)
	nou=lapply(x,function(x) x[x!="u"])
	df=expand.grid(x)
	df.nou=expand.grid(nou)
	return(list(mvms=x,nd=nd,df=df,df.states=df.nou,uncertain=uncertain))
}
