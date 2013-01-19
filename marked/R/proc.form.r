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
		fix.model <- paste("~ ",paste(fix.var, collapse=" + "))
	re.lst <- tms.lst[sapply(tms.lst, "length")==2]
	if(length(re.lst)==0) re.model <- NULL
	else{
		re.model <- lapply(re.lst, function(x){list(model=paste("~",x[1]), sub=paste("~",x[2],"-1", collapse=""))})
	}
	return(list(fix.model=fix.model, re.model=re.model))
}


