#' HMM Observation Probability matrix functions
#' 
#' Functions that compute the probability matrix of the observations given the state for various models. 
#' Currently only CJS, MS models and MS models with state uncertainty are included.
#'  
#' @param pars list of real parameter matrices (id by occasion) for each type of parameter
#' @param m number of states
#' @param T number of occasions
#' @usage cjs_dmat(pars,m,T)
#'        ms_dmat(pars,m,T)
#'        ums_dmat(pars,m,T)
#' @aliases cjs_dmat ms_dmat ums_dmat
#' @return 4-d array of id and occasion-specific observation probability matrices - state-dependent distributions in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_dmat=function(pars,m,T) 
{
	# create 4-d array with a matrix for each id and occasion
	# from pars$p which is a matrix of id by occasion capture probabilities 
	aaply(pars$p,c(1,2),function(p) 
				matrix(c(1-p,1,p,0),nrow=2,ncol=2,byrow=TRUE))
}
ms_dmat=function(pars,m,T) 
{
	# create 4-d array with a matrix for each id and occasion from
	# from pars$p which is a matrix of id by occasion x state capture probabilities 
	# which is split across occasions for multiple states; 
	# each dmat has m+1 rows (0 + m states) and m+1 columns - m states + dead
	aaply(pars$p,1,function(x) 
			{
				laply(split(x,rep(1:(T-1),each=(m-1))),function(p) 
						{
							pdiag=diag(p)
							cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))
						})
			})
}
ums_dmat=function(pars,m,T) 
{
	# create 4-d array with a p matrix for each id and occasion 
	# from pars$p which is a matrix of id by occasion x state capture probabilities 
	# which is split across occasions for multiple states; 
	# also and a 4-d array with delta and
	# then return their product; splits across occasion for multiple states
	pmat=aaply(pars$p,1,function(x) {laply(split(x,rep(1:(T-1),each=(m-1))),function(p) {
							pdiag=diag(p)
							cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))
						})})
	# create 4-d array with a delta matrix for each id and occasion 
	# from pars$delta which is a matrix of id by occasion x state probabilities of identifying state 
	# which is split across occasions for multiple states; 
	deltamat=aaply(pars$delta,1,function(x) {laply(split(x,rep(1:(T-1),each=(m-1))),function(delta) {
				deltamat=matrix(1,ncol=length(delta),nrow=length(delta))
				diag(deltamat)=delta
				return(cbind(rbind(rep(1,ncol(deltamat)),deltamat,1-delta),rep(1,length(delta)+2)))
			})})
    # return the product of the 4-d arrays
    return(pmat*deltamat)
}
