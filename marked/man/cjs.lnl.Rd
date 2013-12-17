\name{cjs.lnl}
\alias{cjs.lnl}
\title{Likelihood function for Cormack-Jolly-Seber model}
\usage{
  cjs.lnl(par, model_data, Phi.links = NULL,
    p.links = NULL, debug = FALSE, all = FALSE, cjsenv)
}
\arguments{
  \item{par}{vector of parameter values}

  \item{model_data}{a list that contains: 1)imat-list of
  vectors and matrices constructed by
  \code{\link{process.ch}} from the capture history data,
  2)Phi.dm design matrix for Phi constructed by
  \code{\link{create.dm}}, 3)p.dm design matrix for p
  constructed by \code{\link{create.dm}}, 4)Phi.fixed
  matrix with 3 columns: ch number(i), occasion number(j),
  fixed value(f) to fix phi(i,j)=f, 5)p.fixed matrix with 3
  columns: ch number(i), occasion number(j), and 6)
  time.intervals intervals of time between occasions if not
  all 1 fixed value(f) to fix p(i,j)=f}

  \item{Phi.links}{vector of links for each parameter}

  \item{p.links}{vector of links for each parameter}

  \item{debug}{if TRUE will printout values of \code{par}
  and function value}

  \item{all}{if TRUE, returns entire list rather than just
  lnl; can be used to extract reals}

  \item{cjsenv}{environment for cjs to update iteration
  counter}
}
\value{
  either -log likelihood value if \code{all=FALSE} or the
  entire list contents of the call to the FORTRAN
  subroutine if \code{all=TRUE}. The latter is used from
  \code{\link{cjs}} after optimization to extract the real
  parameter estimates at the final beta values.
}
\description{
  For a given set of parameters and data, it computes -log
  Likelihood value.
}
\details{
  This function uses a FORTRAN subroutine (cjs.f) to speed
  up computation of the likelihood but the result can also
  be obtained wholly in with a small loss in precision. See
  example for R code. The R and FORTRAN code uses the
  likelihood formulation of Pledger et al.(2003).
}
\examples{
get.p=function(beta,dm,nocc,Fplus)
{
# compute p matrix from parameters (beta) and list of design matrices (dm)
# created by function create.dm
  ps=cbind(rep(1,nrow(dm)/(nocc-1)),
     matrix(dm\%*\%beta,ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE))
  ps[Fplus==1]=plogis(ps[Fplus==1])
  return(ps)
}
get.Phi=function(beta,dm,nocc,Fplus)
{
# compute Phi matrix from parameters (beta) and list of design matrices (dm)
# created by function create.dm
  Phis=cbind(rep(1,nrow(dm)/(nocc-1)),
     matrix(dm\%*\%beta,ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE))
  Phis[Fplus==1]=plogis(Phis[Fplus==1])
  return(Phis)
}
#################################################################################
# cjs.lnl - computes likelihood for CJS using Pledger et al (2003)
# formulation for the likelihood. This code does not cope with fixed parameters or
# loss on capture but could be modified to do so. Also, to work directly with cjs.r and
# cjs.accumulate call to process.ch would have to have all=TRUE to get Fplus and L.
# Arguments:
# par             - vector of beta parameters
# imat           - list of freq, indicator vector and matrices for ch data created by process.ch
# Phi.dm         - list of design matrices; a dm for each capture history
# p.dm           - list of design matrices; a dm for each capture history
# debug          - if TRUE show iterations with par and -2lnl
# time.intervals - intervals of time between occasions
# Value: -LnL using
#################################################################################
cjs.lnl=function(par,model_data,Phi.links=NULL,p.links=NULL,debug=FALSE,all=FALSE) {
	if(debug)cat("\\npar = ",par)
	#extract Phi and p parameters from par vector
	nphi=ncol(model_data$Phi.dm)
	np=ncol(model_data$p.dm)
	beta.phi=par[1:nphi]
	beta.p=par[(nphi+1):(nphi+np)]
	#construct parameter matrices (1 row for each capture history and a column
 #for each occasion)
	Phis=get.Phi(beta.phi,model_data$Phi.dm,nocc=ncol(model_data$imat$chmat),
                          model_data$imat$Fplus)
	if(!is.null(model_data$time.intervals))
	{
		exponent=cbind(rep(1,nrow(Phis)),model_data$time.intervals)
		Phis=Phis^exponent
	}
	ps=get.p(beta.p,model_data$p.dm,nocc=ncol(model_data$imat$chmat),
           model_data$imat$Fplus)
	if(debug)cat("\\npar = ",par)
	# Compute probability of dying in interval from Phis
	M=cbind((1-Phis)[,-1],rep(1,nrow(Phis)))
	# compute cummulative survival from release across each subsequent time
	# and the cummulative probability for detection (capture) across each time
	Phi.cumprod=1-model_data$imat$Fplus + Phis*model_data$imat$Fplus
	cump=(1-model_data$imat$Fplus)+model_data$imat$Fplus*
         (model_data$imat$chmat*ps+(1-model_data$imat$chmat)*(1-ps))
	for (i in 2:ncol(cump))
	{
		Phi.cumprod[,i]=Phi.cumprod[,i-1]*Phi.cumprod[,i]
		cump[,i]=cump[,i-1]*cump[,i]
	}
	# compute prob of capture-history
	pch=rowSums(model_data$imat$L*M*Phi.cumprod*cump)
	lnl=-sum(model_data$imat$freq*log(pch))
	if(debug)cat("\\n-2lnl = ",2*lnl)
	return(lnl)
}
}
\author{
  Jeff Laake
}
\references{
  Pledger, S., K. H. Pollock, et al. (2003). Open
  capture-recapture models with heterogeneity: I.
  Cormack-Jolly-Seber model. Biometrics 59(4):786-794.
}

