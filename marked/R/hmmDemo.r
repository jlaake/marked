#' HMM computation demo functions
#' 
#' Uses fitted hmm model to construct HMM state vectors alpha and phi for demonstration purposes
#' @param model fitted hmm model
#' @param state.names names for states used to label output; if NULL uses strata.labels + Dead state
#' @param obs.names names for observations used to label output; if NULL uses ObsLevels
#' @return hmm demo list which includes 1) lnl - the log-likelihood value, 2) alpha - forward probabilities,
#' 3) beta - backward probabilities, 4) phi - scaled forward probabilities, 5) v- intermediate calculation for phi,
#' 6) dmat - 3-d array with observation probability matrix for each occasion, 7) gamma - 3-d array with state transition probability
#' matrix for each occasion, 8) stateprob - predicted state probabilities, 9) chfowardstrings - set of capture history forward strings, 10) chbackwardstrings - same thing backwards.
#' @author Jeff Laake
#' @export hmmDemo
#' @useDynLib marked
#' @keywords models
#' @examples
#' 
#' # cormack-jolly-seber model
#' data(dipper)
#' mod=crm(dipper,model="hmmcjs")
#' x=hmmDemo(mod,state.names=c("Alive","Dead"),obs.names=c("Missed","Seen"))
#' par(mfrow=c(2,1))
#' barplot(t(x$alpha[45,]),beside=TRUE,names.arg=x$chforwardstrings)
#' barplot(t(x$phi[45,]),beside=TRUE,names.arg=x$chforwardstrings)
#' # multi-state example showing state predictions
#' data(mstrata,package="RMark")
#' mod=crm(mstrata,model="hmmMSCJS")
#' # state predictions are normalized by likelihood value which = rowSums(alpha*beta)
#' cat(paste("\nrowsums = ",rowSums(x$alpha[45,]*x$beta[45,])[1],
#'    "which matches likelihood value",exp(x$lnl[45]),"\n"))
#' # state predictions given the data
#' x$stateprob 
hmmDemo <- function(object,state.names=NULL,obs.names=NULL,...)
{
	result=loglikelihood(object)
	result$beta=backward_prob(object)
	result$stateprob=result$alpha*result$beta/exp(result$lnl)
	if(is.null(state.names))state.names=c(object$data$strata.labels,"Dead")
	if(is.null(obs.names))obs.names=object$data$ObsLevels
    dimnames(result$alpha)[2]=list(state.names)
	dimnames(result$phi)=dimnames(result$alpha)
	dimnames(result$v)=dimnames(result$alpha)
	dimnames(result$beta)=dimnames(result$alpha)
	dimnames(result$stateprob)=dimnames(result$alpha)
	dimnames(result$gamma)[1:2]=list(state.names,state.names)
	names(dimnames(result$gamma))=c("From_state","To_state","Occasion")
	dimnames(result$dmat)[1:2]=list(obs.names,state.names)
	names(dimnames(result$dmat))=c("Observation","State","Occasion")
	colnames(result$alpha)=state.names
	colnames(result$phi)=state.names
	names(dimnames(result$alpha))=c("Occasion","State")
	names(dimnames(result$phi))=names(dimnames(result$alpha))
	names(dimnames(result$beta))=names(dimnames(result$alpha))
	names(dimnames(result$stateprob))=names(dimnames(result$alpha))
	names(dimnames(result$v))=names(dimnames(result$alpha))
	ch=strsplit(object$data$ch[object$data$id==id],",")[[1]]
	result$chfowardstrings=sapply(1:length(ch),function(x)paste(ch[1:x],collapse=""))
	result$chbackwardstrings=sapply(1:length(ch),function(x)paste(rev(ch)[1:x],collapse=""))
	return(result)
}
