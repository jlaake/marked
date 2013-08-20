#' HMM Observation Probability matrix functions
#' 
#' Functions that compute the probability matrix of the observations given the state for various models. 
#' Currently only CJS, MS models and MS models with state uncertainty are included.
#'  
#' @param pars list of real parameter matrices (id by occasion) for each type of parameter
#' @param m number of states
#' @param T number of occasions
#' @aliases cjs_dmat ms_dmat ums_dmat
#' @return 4-d array of id and occasion-specific observation probability matrices - state-dependent distributions in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_dmat=function(pars,m,T) 
{
	# add first occasion p=1
	pmat=array(NA,c(nrow(pars$p),T,2,2))
	for (i in 1:nrow(pmat))
	{
		pmat[i,1,,]=matrix(c(0,1,1,0),nrow=2,ncol=2,byrow=TRUE)
        for(j in 1:(T-1))
        {
			p=pars$p[i,j]
			pmat[i,j+1,,]=matrix(c(1-p,1,p,0),nrow=2,ncol=2,byrow=TRUE)
		}
	}
	pmat
}
ms_dmat=function(pars,m,T) 
{
	# create 4-d array with a matrix for each id and occasion from
	# from pars$p which is a matrix of id by occasion x state capture probabilities 
	# which is split across occasions for multiple states; 
	# each dmat has m+1 rows (0 + m states) and m+1 columns - m states + dead
	# add first occasion p=1
	pmat=array(NA,c(nrow(pars$p),T,m,m))
	for (i in 1:nrow(pmat))
	{
		pdiag=diag(rep(1,m-1))
		pmat[i,1,,]=cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))		
		for(j in 1:(T-1))
		{
			pdiag=diag(pars$p[i,((j-1)*(m-1)+1):(j*(m-1))])			
			pmat[i,j+1,,]=cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))	
		}
	}
	pmat
}
ums_dmat=function(pars,m,T) 
{
	# create 4-d array with a p matrix for each id and occasion 
	# from pars$p which is a matrix of id by occasion x state capture probabilities 
	# which is split across occasions for multiple states; 
	# also and a 4-d array with delta and
	# then return their product; splits across occasion for multiple states
	# add first occasion p=1
	pmat=array(NA,c(nrow(pars$p),T,m+1,m))
	for (i in 1:nrow(pmat))
	{
		pdiag=diag(rep(1,m-1))
		pmat[i,1,,]=cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))		
		for(j in 1:(T-1))
		{
			pdiag=diag(pars$p[i,((j-1)*(m-1)+1):(j*(m-1))])			
			pmat[i,j+1,,]=cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))	
		}
	}
	# create 4-d array with a delta matrix for each id and occasion 
	# from pars$delta which is a matrix of id by occasion-state of 
	# probabilities of identifying state
	# which is split across occasions for multiple states;
	# add first occasion delta=1 for known state at release
	 deltamat=array(NA,c(nrow(pars$delta),T,m+1,m))
	 for (i in 1:nrow(deltamat))
	 {
		 delta=rep(1,m-1)
		 deltax=matrix(1,ncol=length(delta),nrow=length(delta))
		 diag(deltax)=delta
		 deltamat[i,1,,]=cbind(rbind(rep(1,ncol(deltax)),deltax,1-delta),rep(1,length(delta)+2))
		 for(j in 1:(T-1))
		 {
			 delta=pars$delta[i,((j-1)*(m-1)+1):(j*(m-1))]			
			 deltax=matrix(1,ncol=length(delta),nrow=length(delta))
			 diag(deltax)=delta
			 deltamat[i,j+1,,]=cbind(rbind(rep(1,ncol(deltax)),deltax,1-delta),rep(1,length(delta)+2))
		 }
	 }
	 # return the product of the 4-d arrays
    return(pmat*deltamat)
}
ums2_dmat=function(pars,m,T) 
{
	# create 4-d array with a p matrix for each id and occasion 
	# from pars$p which is a matrix of id by occasion x state capture probabilities 
	# which is split across occasions for multiple states; 
	# also and a 4-d array with delta and
	# then return their product; splits across occasion for multiple states
	# add first occasion p=1
	pmat=array(NA,c(nrow(pars$p),T,m+1,m))
	for (i in 1:nrow(pmat))
	{
		pdiag=diag(rep(1,m-1))
		pmat[i,1,,]=cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))		
		for(j in 1:(T-1))
		{
			pdiag=diag(pars$p[i,((j-1)*(m-1)+1):(j*(m-1))])			
			pmat[i,j+1,,]=cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))	
		}
	}
	# create 4-d array with a delta matrix for each id and occasion 
	# from pars$delta which is a matrix of id by occasion-state of 
	# probabilities of identifying state
	# which is split across occasions for multiple states;
	# add first occasion delta=1 for known state at release
	deltamat=array(NA,c(nrow(pars$delta),T,m+1,m))
	for (i in 1:nrow(deltamat))
	{
		delta=rep(1,m-1)
		deltax=matrix(1,ncol=length(delta),nrow=length(delta))
		diag(deltax)=delta
		deltamat[i,1,,]=cbind(rbind(rep(1,ncol(deltax)),deltax,1-delta),rep(1,length(delta)+2))
		for(j in 1:(T-1))
		{
			delta=pars$delta[i,((j-1)*(m-1)+1):(j*(m-1))]			
			deltax=matrix(1,ncol=length(delta),nrow=length(delta))
			diag(deltax)=delta
			deltamat[i,j+1,,]=cbind(rbind(rep(1,ncol(deltax)),deltax,1-delta),rep(1,length(delta)+2))
		}
	}
	# return the product of the 4-d arrays
	return(pmat*deltamat)
}
