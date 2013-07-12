#' HMM Observation Probability matrix functions
#' 
#' Functions that compute the probability matrix of the observations given the state for various models. 
#' Currently only CJS, MS models and MS models with state uncertainty are included.
#'  
#' @param id sequential id for the observed sequence
#' @param ddl design data list of parameters for each id
#' @param parlist list of parameter vectors split by type (eg Phi, p in CJS)
#' @param parameters formulas for each parameter type
#' @usage cjs_dmat(pars,m,T)
#'        ms_dmat(pars,m,T)
#'        ums_dmat(pars,m,T)
#' @aliases cjs_dmat ms_dmat ums_dmat
#' @return array of id and occasion-specific observation probability matrices - state-dependent distributions in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_dmat=function(pars,m,T) 
{
	# Arguments:
	# parlist: list of parameter vectors split by type (eg Phi, p in CJS)
	# parameters: formulas for each parameter type
	return(aaply(pars$p,c(1,2),function(p) matrix(c(1-p,1,p,0),nrow=2,ncol=2,byrow=TRUE)))
}
ms_dmat=function(pars,m,T) 
{
	# create array with a matrix for each id and occasion
	return(aaply(pars$p,1,function(x) {laply(split(x,rep(1:(T-1),each=(m-1))),function(p) {
						pdiag=diag(p)
						cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))
					})}))
}
ums_dmat=function(pars,m,T) 
{
	# create array with a p matrix for each id and occasion and one with delta and
	# then return their product
	pmat=aaply(pars$p,1,function(x) {laply(split(x,rep(1:(T-1),each=(m-1))),function(p) {
							pdiag=diag(p)
							cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))
						})})
	deltamat=aaply(pars$delta,1,function(x) {laply(split(x,rep(1:(T-1),each=(m-1))),function(delta) {
				deltamat=matrix(1,ncol=length(delta),nrow=length(delta))
				diag(deltamat)=delta
				return(cbind(rbind(rep(1,ncol(deltamat)),deltamat,1-delta),rep(1,length(delta)+2)))
			})})
	return(pmat*deltamat)
}
