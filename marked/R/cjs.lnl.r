#' Likelihood function for Cormack-Jolly-Seber model
#' 
#' For a given set of parameters and data, it computes -2*log Likelihood value.
#' 
#' This function uses a FORTRAN subroutine to compute the likelihood but the
#' result can also be obtained from the following functions written wholly in
#' R.  The code uses the likelihood formulation of Pledger et al.(2003).
#' 
#' \preformatted{ # compute p matrix from parameters (beta) and list of design
#' matrices (dm) # created by function create.dm
#' get.p=function(beta,dm,nocc,Fplus) {
#' ps=cbind(rep(1,nrow(dm)/(nocc-1)),matrix(dm%*%beta,ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE))
#' ps[Fplus==1]=plogis(ps[Fplus==1]) return(ps) }
#' 
#' # compute Phi matrix from parameters (beta) and list of design matrices (dm)
#' # created by function create.dm get.Phi=function(beta,dm,nocc,Fplus) {
#' Phis=cbind(rep(1,nrow(dm)/(nocc-1)),matrix(dm%*%beta,ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE))
#' Phis[Fplus==1]=plogis(Phis[Fplus==1]) return(Phis) }
#' 
#' #################################################################################
#' # cjs.lnl - computes likelihood for CJS using Pledger et al (2003)
#' formulation # does not cope with fixed parameters or loss on capture but
#' could be # modified to do so. # # Arguments: # # par - vector of beta
#' parameters # imat - list of freq, indicator vector and matrices for ch data
#' created by process.ch # Phi.dm - list of design matrices; a dm for each
#' capture history # p.dm - list of design matrices; a dm for each capture
#' history # debug - if TRUE show iterations with par and -2lnl #
#' time.intervals - intervals of time between occasions # # Value: -2LnL using
#' Pledger et al (2003) formulation for the likelihood
#' #################################################################################
#' cjs.lnl=function(par,imat,Phi.dm,p.dm,debug=FALSE,time.intervals=NULL) {
#' if(debug)cat("\npar = ",par)
#' ###################################################################################
#' # compute Phi and p matrices (a value for each occasion for each unique ch
#' (or animal))
#' ###################################################################################
#' #extract Phi and p parameters from par vector nphi=ncol(Phi.dm)
#' np=ncol(p.dm) beta.phi=par[1:nphi] beta.p=par[(nphi+1):(nphi+np)] #
#' construct parameter matrices (1 row for each capture history and a column
#' for each occasion)
#' Phis=get.Phi(beta.phi,Phi.dm,nocc=ncol(imat$chmat),imat$Fplus)
#' if(!is.null(time.intervals)) {
#' exponent=matrix(c(1,time.intervals),nrow=nrow(Phis),ncol=ncol(Phis),byrow=TRUE)
#' Phis=Phis^exponent } ps=get.p(beta.p,p.dm,nocc=ncol(imat$chmat),imat$Fplus)
#' ###################################################################################
#' # construct cumphi and cump for known portion of ch from release to last
#' recapture
#' ###################################################################################
#' # compute cummulative survival from release across each subsequent time
#' Phi.cumprod=cumprod.matrix(1-imat$Fplus + Phis*imat$Fplus) # extract the
#' cummulative Phi for each individual (value on last occasion)
#' cumPhi=Phi.cumprod[cbind(1:nrow(imat$chmat),imat$last)] # compute the
#' cummulative probability for detection (capture)
#' cump=prod.matrix((1-imat$FtoL)+imat$FtoL*(imat$chmat*ps+(1-imat$chmat)*(1-ps)))
#' ###################################################################################
#' # construct likelihood portion for 0...0 tail (chi in typical formulation)
#' ###################################################################################
#' # first construct a matrix with potential mortality in each column after
#' last occasion seen # the last column is all 1's for case where animal
#' survived to end of study
#' mort.last=cbind((1-imat$Lplus+imat$Lplus*(1-Phis))[,-1],rep(1,nrow(Phis))) #
#' next construct the cummulative product of 1-p for each occasion alive and
#' not seen in tail given alive p.last=cumprod.matrix(1-imat$Lplus +
#' imat$Lplus*(1-ps))*imat$L # next compute the sum over all possible times of
#' death after last occasion seen; this includes the case where # the animal is
#' still alive at last occasion p0=rowSums(Phi.cumprod/cumPhi*p.last*mort.last)
#' ###################################################################################
#' # return -2LnL which is the sum of the logarithm of the 3 components
#' ###################################################################################
#' lnl=-2*(sum(imat$freq*log(p0)) + sum(imat$freq*log(cump)) +
#' sum(imat$freq*log(cumPhi))) if(debug)cat("\nlnl = ",lnl) return(lnl) } }
#' 
#' @param par vector of parameter values
#' @param model_data a list that contains: 1)imat-list of vectors and matrices constructed by
#' \code{\link{process.ch}} from the capture history data, 2)Phi.dm design matrix for Phi constructed by \code{\link{create.dm}},
#' 3)p.dm design matrix for p constructed by \code{\link{create.dm}},
#' 4)Phi.fixed matrix with 3 columns: ch number(i), occasion number(j),
#' fixed value(f) to fix phi(i,j)=f, 5)p.fixed matrix with 3 columns: ch number(i), occasion number(j), and
#' 6) time.intervals intervals of time between occasions if not all 1
#' fixed value(f) to fix p(i,j)=f
#' @param Phi.links vector of links for each parameter
#' @param p.links vector of links for each parameter
#' @param debug if TRUE will printout values of \code{par} and function value
#' @param all if TRUE, returns entire list rather than just lnl; can be used to
#' extract reals
#' @return either -2*log likelihood value if \code{all=FALSE} or the entire
#' list contents of the call to the FORTRAN subroutine if \code{all=TRUE}. The
#' latter is used from \code{\link{cjs}} after optimization to extract the real
#' parameter estimates at the final beta values.
#' @author Jeff Laake
#' @references Pledger, S., K. H. Pollock, et al. (2003). Open
#' capture-recapture models with heterogeneity: I. Cormack-Jolly-Seber model.
#' Biometrics 59(4):786-794.
cjs.lnl=function(par,model_data,Phi.links=NULL,p.links=NULL,debug=FALSE,all=FALSE)
{
	f_eval=get(".markedfunc_eval",envir=.GlobalEnv)+1
	assign(".markedfunc_eval", f_eval, envir = .GlobalEnv)
	if(debug)cat("\npar = ",par)
	nocc=model_data$imat$nocc
	nphi=ncol(model_data$Phi.dm)
	np=ncol(model_data$p.dm)
	beta.phi=par[1:nphi]
	beta.p=par[(nphi+1):(nphi+np)]
	Phibeta=as.vector(model_data$Phi.dm%*%beta.phi)
#	if(length(Phi.links>0))
#	{
#		Phibeta[Phi.links]=as.vector((sin(Phi.dm[Phi.links,,drop=FALSE]%*%beta.phi)+1)/2)
#	    Phibeta[Phi.links]=log(1/(1/Phibeta[Phi.links]-1))
#	}
	Phibeta=matrix(Phibeta,ncol=nocc-1,nrow=nrow(model_data$Phi.dm)/(nocc-1),byrow=TRUE)
	pbeta=as.vector(model_data$p.dm%*%beta.p)
#	if(length(p.links>0))
#	{
#	   pbeta[p.links]=as.vector((sin(p.dm[p.links,,drop=FALSE]%*%beta.p)+1)/2)
#	   pbeta[p.links]=log(1/(1/pbeta[p.links]-1))
#    }
	pbeta=matrix(pbeta,ncol=nocc-1,nrow=nrow(model_data$p.dm)/(nocc-1),byrow=TRUE)	
	lnl=.Fortran("cjs",as.double(model_data$imat$chmat),as.double(Phibeta),as.double(pbeta),
			as.double(model_data$imat$first),as.double(model_data$imat$last),as.double(model_data$imat$freq),
			as.integer(model_data$imat$loc),as.double(model_data$Phi.fixed),as.double(model_data$p.fixed),
			as.double(model_data$time.intervals),as.integer(nrow(model_data$imat$chmat)),           
			as.integer(ncol(model_data$imat$chmat)),as.integer(nrow(model_data$Phi.fixed)),
			as.integer(nrow(model_data$p.fixed)),lnl=double(1),p0=double(nrow(model_data$imat$chmat)),PACKAGE="marked")
	if(debug)
	{
		cat("\n-2lnl = ",lnl$lnl)
	} else
	if((f_eval-100*floor(f_eval/100))==0)
	{
	    cat("\r Number of evaluations: ",f_eval," -2lnl:",formatC(lnl$lnl,digits=10))
		flush.console()
	}	
	if(all)
		return(lnl)
	else
		return(lnl$lnl)
}
