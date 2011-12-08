#' Likelihood function for Jolly-Seber model using Schwarz-Arnason POPAN
#' formulation
#' 
#' For a given set of parameters and data, it computes -2*log Likelihood value
#' but does not include data factorials. Factorials for unmarked are not needed
#' but are included in final result by \code{\link{js}} so the result matches
#' output from MARK for the POPAN model.
#' 
#' This functions uses \code{\link{cjs.lnl}} and then supplements with the
#' remaining calculations to compute the likelihood for the POPAN formulation
#' (Arnason and Schwarz 1996) of the Jolly-Seber model.
#' 
#' @param par vector of parameter values
#' @param imat list of vectors and matrices constructed by
#' \code{\link{process.ch}} from the capture history data.
#' @param Phi.dm design matrix for Phi constructed by \code{\link{create.dm}}
#' @param p.dm design matrix for p constructed by \code{\link{create.dm}}
#' @param pent.dm design matrix for probability of entry constructed by
#' \code{\link{create.dm}}
#' @param N.dm design matrix for estimates of number of animals not caught from
#' super-population constructed by \code{\link{create.dm}}
#' @param Phi.fixed matrix with 3 columns: ch number(i), occasion number(j),
#' fixed value(f) to fix phi(i,j)=f
#' @param p.fixed matrix with 3 columns: ch number(i), occasion number(j),
#' fixed value(f) to fix p(i,j)=f
#' @param pent.fixed matrix with 3 columns: ch number(i), occasion number(j),
#' fixed value(f) to fix pent(i,j)=f
#' @param debug if TRUE will printout values of \code{par} and function value
#' @param time.intervals intervals of time between occasions if not all 1
#' @param nobstot number of unique caught at least once by group if applicable
#' @return -2*log likelihood value (excluding data factorials)
#' @author Jeff Laake
#' @references Schwarz, C. J., and A. N. Arnason. 1996. A general methodology
#' for the analysis of capture-recapture experiments in open populations.
#' Biometrics 52:860-873.
js.lnl=function(par,imat,Phi.dm,p.dm,pent.dm,N.dm,Phi.fixed=NULL,p.fixed=NULL,
                 pent.fixed=NULL,debug=FALSE,time.intervals=NULL,nobstot)
###################################################################################
# js.lnl - computes -2*log likelihood value of js model using call to cjs.lnl 
#            and adds on the entry portion of the js likelihood
#
# Arguments:
#    par             - vector of parameter values to evaluate likelihood
#    imat            - list of indicator vectors, frequency etc from process.ch
#    Phi.dm          - Phi design matrix created by call to create.dm
#    p.dm            - p design matrix created by call to create.dm
#    pent.dm         - probability of entry desgn matrix
#    N.dm            - super-population size design matrix
#    Phi.fixed       - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix phi(i,j)=f 
#    p.fixed         - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
#    pent.fixed      - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
#    debug           - if TRUE show iterations with par and -2lnl
#    time.intervals  - intervals of time between occasions
#    nobstot         - number of unique caught at least once by group if applicable
#
# Value: -2*log likelihood
###################################################################################
{
	get.pent=function(beta,dm,nocc)
	{
		pents=cbind(rep(1,nrow(dm)/(nocc-1)),exp(matrix(as.vector(dm%*%beta),ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE)))
		pents=pents/apply(pents,1,sum)
		return(pents)
	}
	get.p=function(beta,dm,nocc)
	{
		ps=cbind(rep(1,nrow(dm)/(nocc-1)),plogis(matrix(as.vector(dm%*%beta),ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE)))
		return(ps)
	}
# compute Phi matrix from parameters (beta) and list of design matrices (dm)
# created by function create.dm
	get.Phi=function(beta,dm,nocc)
	{
		Phis=cbind(rep(1,nrow(dm)/(nocc-1)),plogis(matrix(as.vector(dm%*%beta),ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE)))
		return(Phis)
	}
	
if(debug)cat("\npar = ",par)
f_eval=get(".markedfunc_eval",envir=.GlobalEnv)+1
assign(".markedfunc_eval", f_eval, envir = .GlobalEnv)
# initialize constants and parameter vectors
nocc=imat$nocc
nphi=ncol(Phi.dm)
np=ncol(p.dm)
npent=ncol(pent.dm)
nN=ncol(N.dm)
beta.phi=par[1:nphi]
beta.p=par[(nphi+1):(nphi+np)]
beta.pent=par[(nphi+np+1):(nphi+np+npent)]
beta.N=par[(nphi+np+npent+1):(nphi+np+npent+nN)]
# create Phi and p beta matrices excluding first occasion on p
Phibeta=matrix(as.vector(Phi.dm%*%beta.phi),ncol=nocc-1,nrow=nrow(Phi.dm)/(nocc-1),byrow=TRUE)
pbeta=matrix(as.vector(p.dm%*%beta.p),ncol=nocc,nrow=nrow(p.dm)/(nocc),byrow=TRUE)
if(is.null(time.intervals))time.intervals=rep(1,nocc-1)
# compute CJS portion of the likelihood
cjslnl=.Fortran("cjs",as.double(imat$chmat),as.double(Phibeta),as.double(pbeta[,-1]),
           as.double(imat$first),as.double(imat$last),as.double(imat$freq),
           as.integer(imat$loc),as.double(Phi.fixed),as.double(p.fixed[p.fixed[,2]!=1,,drop=FALSE]),
           as.double(time.intervals),as.integer(nrow(imat$chmat)),           
           as.integer(ncol(imat$chmat)),as.integer(nrow(Phi.fixed)),
           as.integer(nrow(p.fixed[p.fixed[,2]!=1,,drop=FALSE])),lnl=double(1),
           p0=double(nrow(imat$chmat)),PACKAGE="marked")
# next add on likelihood component for first capture
pents=get.pent(beta.pent,pent.dm,nocc)
pents.dummy=pents[imat$freq==0,]
ps=plogis(matrix(as.vector(p.dm%*%beta.p),ncol=nocc,nrow=nrow(p.dm)/(nocc),byrow=TRUE))
ps.dummy=ps[imat$freq==0,]
Phis=get.Phi(beta.phi,Phi.dm,nocc)
p.occ=ps[cbind(1:nrow(pents),imat$first)]
ps[,nocc]=0
Phis=cbind(Phis[,2:nocc],rep(1,nrow(Phis)))
entry.p=(1-ps)*Phis*(1-imat$First)+imat$First
entry.p=(t(apply(entry.p,1,function(x) rev(cumprod(rev(x)))))*pents*(1-imat$Fplus))
entry.p=apply(entry.p,1,sum)*p.occ
lnl=cjslnl$lnl-2*sum(imat$freq*log(entry.p))
# next add on likelihood component for those not caught from dummy 1000000,0100000,...data 
# return complete likelihood value except that calling function js adds the ui factorials to match
# POPAN output from MARK
Ns=exp(as.vector(N.dm%*%beta.N))
for (i in 1:length(Ns))
{
  index0=(i-1)*nocc+1
  index1=i*nocc
  ps=ps.dummy[index0:index1,]
  pents=pents.dummy[index0:index1,]
  lnl=lnl-2*sum(Ns[i]*log(sum(diag(pents)*cjslnl$p0[imat$freq==0][index0:index1]*(1-diag(ps)))))
  lnl=lnl-2*lfactorial(nobstot[i]+Ns[i])+2*lfactorial(Ns[i])
}
if(debug)
{
	cat("\nlnl = ",lnl)
} else
    if((f_eval-100*floor(f_eval/100))==0)
    {
		cat("\r Number of evaluations: ",f_eval," -2lnl:",formatC(lnl,digits=10))
		flush.console()
    }	
return(lnl)
}
