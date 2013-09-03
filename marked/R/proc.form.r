#' Mixed effect model formula parser
#' 
#' Parses a mixed effect model in the lme4 structure of  ~fixed +(re1|g1) +...+(ren|gn)
#' 
#' @param f formula for mixed effect mode in the form used in lme4; ~fixed +(re1|g1) +...+(ren|gn)
#' @return A list with elements fix.model and re.model. fix.model contains the formula for the fixed effects;
#' re.model contains elements sub, the grouping formula and model the design formula for the
#' random effect. Each formula is of type character and must be wrapped with as.formula in use with model.matrix
#' @author Devin Johnson <devin.johnson@@noaa.gov>
#' 
proc.form <- function(f){
	tms <- terms(f)
	tms.lab <- attr(tms, "term.labels")
	tms.lst <- strsplit(tms.lab, rep(" | ",length(tms.lab)), fixed=TRUE)
	fix.var <- attr(tms, "term.labels")[sapply(tms.lst, "length")==1]
	if(length(fix.var)==0) 
		fix.model="~1"
	else
	if(length(grep("-1",f))==0)
		fix.model <- paste("~ ",paste(fix.var, collapse=" + "))
	else
		fix.model <- paste("~ -1 +",paste(fix.var, collapse=" + "))
	re.lst <- tms.lst[sapply(tms.lst, "length")==2]
	if(length(re.lst)==0) re.model <- NULL
	else{
		re.model <- lapply(re.lst, function(x){list(model=paste("~",x[1]), sub=paste("~",x[2],"-1", collapse=""))})
		names(re.model) <- sapply(re.lst, function(x) x[2])
	}
	return(list(fix.model=fix.model, re.model=re.model))
}
#' Create design matrices for random effect models
#' 
#' Takes output from \code{proc.form} and creates a list of design marices for use in MCMC sampler 
#' functions.
#' 
#'  @param mlist a named list created from \code{\link{proc.form}}
#'  @return A list containing design matrices for each random effect in the model processed 
#'  by \code{proc.form}.
#'  @author Devin Johnson <devin.johnson@@noaa.gov>

re.design.mat<-function(mlist, data){
	redm.list=vector("list",length(mlist$re))
	if(is.null(mlist$re)) return(NULL)
	else{
		for(i in 1:length(mlist$re)){
			sub=strsplit(mlist$re.model[[i]]$sub, "~ ")[[1]][2]
			mod=strsplit(mlist$re.model[[i]]$mod, "~ ")[[1]][2]
			if(mod=="1") form=formula(paste(c("~", sub), collapse=""))
			else form=formula(paste(c("~", mod, ":(", sub, ")"), collapse=""))
			redm.list[[i]]=Matrix(model.matrix(form, data))
			names(redm.list)[i]=paste(c(mod, "|",strsplit(sub, " -1")), collapse="")
		}
	}
}

