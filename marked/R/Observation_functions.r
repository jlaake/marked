#' HMM Observation Probability matrix functions
#' 
#' Functions that compute the probability matrix of the observations given the state for various models. 
#' Currently only CJS, MS models and MS models with state uncertainty are included.
#'  
#' @param pars list of real parameter matrices (id by occasion) for each type of parameter
#' @param m number of states
#' @param F initial occasion vector 
#' @param T number of occasions
#' @aliases cjs_dmat ms_dmat ums_dmat
#' @return 4-d array of id and occasion-specific observation probability matrices - state-dependent distributions in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_dmat=function(pars,m,F,T) 
{
	value=.Fortran("cjsp",as.double(pars$p),as.integer(nrow(pars$p)),
			as.integer(F),as.integer(T),pmat=double(nrow(pars$p)*T*4),PACKAGE="marked")
	dim(value$pmat)=c(nrow(pars$p),T,2,2)
	value$pmat
}
ms_dmat=function(pars,m,F,T) 
{
	if(is.list(m))m=m$ns*m$na+1
	value=.Fortran("msp",as.double(pars$p),as.integer(nrow(pars$p)),as.integer(m),
			as.integer(F),as.integer(T),pmat=double(nrow(pars$p)*T*m^2),PACKAGE="marked")
	dim(value$pmat)=c(nrow(pars$p),T,m,m)
	value$pmat
}
ums_dmat=function(pars,m,F,T) 
{
	if(is.list(m))m=m$ns*m$na+1
	value=.Fortran("umsp",as.double(pars$p),as.double(pars$delta),as.integer(nrow(pars$p)),as.integer(m),
			as.integer(F),as.integer(T),pmat=double(nrow(pars$p)*T*(m+1)*m),PACKAGE="marked")
	dim(value$pmat)=c(nrow(pars$p),T,m+1,m)
	value$pmat
}
ums2_dmat=function(pars,m,F,T) 
{
	value=.Fortran("ums2p",as.double(pars$p),as.double(pars$delta),as.integer(nrow(pars$p)),as.integer(m$na),
			as.integer(m$ns),as.integer(F),as.integer(T),
			pmat=double(nrow(pars$p)*T*(m$na*(m$ns+1)+1)*(m$ns*m$na+1)),PACKAGE="marked")
	dim(value$pmat)=c(nrow(pars$p),T,m$na*(m$ns+1)+1,m$ns*m$na+1)
	value$pmat
}
# cjs_dmat=function(pars,m,F,T) 
# {
# 	# add first occasion p=1
# 	pmat=array(NA,c(nrow(pars$p),T,2,2))
# 	for (i in 1:nrow(pmat))
# 	{
# 		pmat[i,F[i],,]=matrix(c(0,1,1,0),nrow=2,ncol=2,byrow=TRUE)
#      for(j in F[i]:(T-1))
#         {
# 			p=pars$p[i,j]
# 			pmat[i,j+1,,]=matrix(c(1-p,1,p,0),nrow=2,ncol=2,byrow=TRUE)
# 		}
# 	}
# 	pmat
# }
#ms_dmat=function(pars,m,F,T) 
#{
#	# create 4-d array with a matrix for each id and occasion from
#	# from pars$p which is a matrix of id by occasion x state capture probabilities 
#	# which is split across occasions for multiple states; 
#	# each dmat has m+1 rows (0 + m states) and m+1 columns - m states + dead
#	# add first occasion p=1
#	if(is.list(m))m=m$ns*m$na+1
#	pmat=array(NA,c(nrow(pars$p),T,m,m))
#	for (i in 1:nrow(pmat))
#	{
#		pdiag=diag(rep(1,m-1))
#		pmat[i,F[i],,]=cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))		
#		for(j in F[i]:(T-1))
#		{
#			pdiag=diag(pars$p[i,((j-1)*(m-1)+1):(j*(m-1))])			
#			pmat[i,j+1,,]=cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))	
#		}
#	}
#	pmat
#}
#ums_dmat=function(pars,m,F,T) 
#{
#	# create 4-d array with a p matrix for each id and occasion 
#	# from pars$p which is a matrix of id by occasion x state capture probabilities 
#	# which is split across occasions for multiple states; 
#	# also and a 4-d array with delta and
#	# then return their product; splits across occasion for multiple states
#	# add first occasion p=1
#	pmat=array(NA,c(nrow(pars$p),T,m+1,m))
#	for (i in 1:nrow(pmat))
#	{
#		pdiag=diag(rep(1,m-1))
#		pmat[i,F[i],,]=cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))		
#		for(j in F[i]:(T-1))
#		{
#			pdiag=diag(pars$p[i,((j-1)*(m-1)+1):(j*(m-1))])			
#			pmat[i,j+1,,]=cbind(rbind(1-colSums(pdiag),pdiag,colSums(pdiag)),c(1,rep(0,nrow(pdiag)+1)))	
#		}
#	}
#	# create 4-d array with a delta matrix for each id and occasion 
#	# from pars$delta which is a matrix of id by occasion-state of 
#	# probabilities of identifying state
#	# which is split across occasions for multiple states;
#	# add first occasion delta=1 for known state at release
#	 deltamat=array(NA,c(nrow(pars$delta),T,m+1,m))
#	 for (i in 1:nrow(deltamat))
#	 {
#		 delta=rep(1,m-1)
#		 deltax=matrix(1,ncol=length(delta),nrow=length(delta))
#		 diag(deltax)=delta
#		 deltamat[i,F[i],,]=cbind(rbind(rep(1,ncol(deltax)),deltax,1-delta),rep(1,length(delta)+2))
#		 for(j in F[i]:(T-1))
#		 {
#			 delta=pars$delta[i,((j-1)*(m-1)+1):(j*(m-1))]			
#			 deltax=matrix(1,ncol=length(delta),nrow=length(delta))
#			 diag(deltax)=delta
#			 deltamat[i,j+1,,]=cbind(rbind(rep(1,ncol(deltax)),deltax,1-delta),rep(1,length(delta)+2))
#		 }
#	 }
#	 # return the product of the 4-d arrays
#   return(pmat*deltamat)
#}
#ums2_dmat=function(pars,m,F,T) 
#{
#	# create 4-d array with a p matrix for each id and occasion 
#	# from pars$p which is a matrix of id by occasion x state capture probabilities 
#	# which is split across occasions for multiple states; 
#	# also and a 4-d array with delta and
#	# then return their product; splits across occasion for multiple states
#	# add first occasion p=1
#	pmat=array(0,c(nrow(pars$p),T,m$na*(m$ns+1)+1,m$ns*m$na+1))
#	# loop over rows (individual cap histories)
#	for (i in 1:nrow(pmat))
#	{
#		# first occasion p=1 
#		for(k in 1:m$na)
#		{
#			cols=(m$ns*(k-1)+1):(m$ns*k)
#			rows=((m$ns+1)*(k-1)+2):((m$ns+1)*k)
#			diag(pmat[i,1,rows,cols])=1
#			pmat[i,F[i],k*(m$ns+1)+1,cols]=1
#			
#		}
#		pmat[i,F[i],1,m$ns*m$na+1]=1 # death state
#		# loop over each remaining occasion using estimates of p at each occasion
#		for(j in F[i]:(T-1))
#		{
#			px=pars$p[i,((j-1)*(m$na*m$ns)+1):(j*(m$na*m$ns))]
#			for(k in 1:m$na)
#			{
#				cols=(m$ns*(k-1)+1):(m$ns*k)
#				rows=((m$ns+1)*(k-1)+2):((m$ns+1)*k)	
#			    pvec=px[(m$ns*(k-1)+1):(k*m$ns)]
#				diag(pmat[i,j+1,rows,cols])=pvec
#				pmat[i,j+1,k*(m$ns+1)+1,cols]=pvec
#				pmat[i,j+1,1,cols]=1-pvec
#			}
#			pmat[i,j+1,1,m$ns*m$na+1]=1 # death state
#		}
#	}
#	# create 4-d array with a delta matrix for each id and occasion 
#	# from pars$delta which is a matrix of id by occasion-state of 
#	# probabilities of identifying state
#	# which is split across occasions for multiple states;
#	# add first occasion delta=1 for known state at release
#	deltamat=array(1,c(nrow(pars$delta),T,m$na*(m$ns+1)+1,m$ns*m$na+1))
#	for (i in 1:nrow(deltamat))
#	{
#		for(k in 1:m$na)
#		{
#			cols=(m$ns*(k-1)+1):(m$ns*k)
#			deltamat[i,F[i],k*(m$ns+1)+1,cols]=0
#		}
#		# loop over each remaining occasion using estimates of p at each occasion
#		for(j in F[i]:(T-1))
#		{
#		px=pars$delta[i,((j-1)*(m$na*m$ns)+1):(j*(m$na*m$ns))]
#		for(k in 1:m$na)
#			{
#				cols=(m$ns*(k-1)+1):(m$ns*k)
#				rows=((m$ns+1)*(k-1)+2):((m$ns+1)*k)	
#				pvec=px[(m$ns*(k-1)+1):(k*m$ns)]
#				diag(deltamat[i,j+1,rows,cols])=pvec
#				deltamat[i,j+1,k*(m$ns+1)+1,cols]=1-pvec
#			}
#		}
#	}
#	# return the product of the 4-d arrays
#	return(pmat*deltamat)
#}

