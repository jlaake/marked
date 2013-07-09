#' HMM Transition matrix functions
#' 
#' Functions that compute the transition matrix for various models. Currently only CJS and MS models
#' are included.
#'  
#' @param id sequential id for the observed sequence
#' @param ddl design data list of parameters for each id
#' @param parlist list of parameter vectors split by type (eg Phi, p in CJS)
#' @param parameters formulas for each parameter type
#' @usage cjs_gamma(id,ddl,parlist,parameters)
#'        ms_gamma(id,ddl,parlist,parameters)
#' @aliases cjs_gamma ms_gamma
#' @return list of occasion-specific transition matrices - Gamma in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_gamma=function(id,ddl,parlist,parameters) 
{
	# Arguments:
	# parlist: list of parameter vectors split by type (eg Phi, p in CJS)
	# id: sequential id for the observed sequence
	# ddl: design data list of parameters for each id
	# parameters: formulas for each parameter type
	# 
	# select ddl for Phi and extract those for this id
	phiddl=ddl[["Phi"]]
	phiddl=phiddl[phiddl$id==id,,drop=FALSE]
	# Using the formula and design data for Phi(survival), create the design matrix
	# with model.matrix, then multiply design matrix and parameters and 
	# use inverse logit function plogis to compute Phi
	phis=reals("Phi",ddl=phiddl,parameters=parameters,link="logit",parlist=parlist)
	# create list with a transition matrix for each occasion for CJS
	return(lapply(phis,function(phi) matrix(c(phi,1-phi,0,1),nrow=2,ncol=2,byrow=T)))
}
ms_gamma=function(id,ddl,parlist,parameters) 
{
	# select ddl for Phi and Psi and extract those for this id 
	phiddl=ddl[["Phi"]]
	phiddl=phiddl[phiddl$id==id,,drop=FALSE]
	psiddl=ddl[["Psi"]]
	psiddl=psiddl[psiddl$id==id,,drop=FALSE]
	# compute values of Phi for each occasion
	phis=reals("Phi",ddl=phiddl,parameters=parameters,link="logit",parlist=parlist)
	# compute values of Psi for each occasion
	psis=reals("Psi",ddl=psiddl,parameters=parameters,link="log",parlist=parlist)
	# create list with a matrix for each occasion
	phimat=tapply(phis,phiddl$time,function(phi) {
				phimat=matrix(phi,ncol=length(phi),nrow=length(phi))
				return(rbind(cbind(phimat,1-phi),c(rep(0,length(phi)),1)))
			})
	psimat=tapply(psis,psiddl$time,function(psi) {
				psimat=matrix(psi,ncol=sqrt(length(psi)),byrow=TRUE)
				psimat=psimat/rowSums(psimat)
				return(rbind(cbind(psimat,rep(1,nrow(psimat))),rep(1,nrow(psimat)+1)))
			})
	return(lapply(1:length(phimat), function(i) psimat[[i]]*phimat[[i]]))
}
