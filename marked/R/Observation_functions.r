#' HMM Observation Probability matrix functions
#' 
#' Functions that compute the probability matrix of the observations given the state for various models. 
#' Currently only CJS, MS models and MS models with state uncertainty are included.
#'  
#' @param id sequential id for the observed sequence
#' @param ddl design data list of parameters for each id
#' @param parlist list of parameter vectors split by type (eg Phi, p in CJS)
#' @param parameters formulas for each parameter type
#' @usage cjs_dmat(id,ddl,parlist,parameters)
#'        ms_dmat(id,ddl,parlist,parameters)
#'        ums_dmat(id,ddl,parlist,parameters)
#' @aliases cjs_dmat ms_dmat ums_dmat
#' @return list of occasion-specific observation probability matrices - state-dependent distributions in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_dmat=function(id,ddl,parlist,parameters) 
{
	# Arguments:
	# parlist: list of parameter vectors split by type (eg Phi, p in CJS)
	# id: sequential id for the observed sequence
	# ddl: design data list of parameters for each id
	# parameters: formulas for each parameter type
	# 
	# select ddl for p and extract those for this id
	pddl=ddl[["p"]]
	pddl=pddl[pddl$id==id,,drop=FALSE]
	# Using the formula and design data for p, create the design matrix
	# with model.matrix, then multiply design matrix and parameters and 
	# use inverse logit function plogis to compute p
	ps=reals("p",ddl=pddl,parameters=parameters,link="logit",parlist=parlist)
	# create list with a matrix for each occasion containing D for CJS
	return(lapply(ps,function(p) matrix(c(1-p,1,p,0),nrow=2,ncol=2,byrow=T)))
}
ms_dmat=function(id,ddl,parlist,parameters) 
{
	# select ddl for p and extract those for this id
	pddl=ddl[["p"]]
	pddl=pddl[pddl$id==id,,drop=FALSE]
	# compute values of p for each occasion
	ps=reals("p",ddl=pddl,parameters=parameters,link="logit",parlist=parlist)
	# create list with a matrix for each occasion
	return(tapply(ps,pddl$time,function(p) {
						pdiag=diag(p)
						cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))
					}))
}
ums_dmat=function(id,ddl,parlist,parameters) 
{
	# Arguments:
	# parlist: list of parameter vectors split by type (eg Phi, p, Psi)
	# id: sequential id for the observed sequence
	# ddl: design data list of parameters for each id
	# parameters: formulas for each parameter type
	# 
	# select ddl for p and delta and extract those for this id
	pddl=ddl[["p"]]
	pddl=pddl[pddl$id==id,,drop=FALSE]
	deltaddl=ddl[["delta"]]
	deltaddl=deltaddl[deltaddl$id==id,,drop=FALSE]
	# compute values of p for each occasion/stratum
	ps=reals("p",ddl=pddl,parameters=parameters,link="logit",parlist=parlist)
	# compute values of delta for each occasion/stratum
	deltas=reals("delta",ddl=deltaddl,parameters=parameters,link="logit",parlist=parlist)
	# create list with a p matrix for each occasion and one with delta and
	# then return their product
	pmat=tapply(ps,pddl$time,
			function(p) 
			{
				pdiag=diag(p)
				cbind(rbind(1-colSums(pdiag),pdiag),c(1,rep(0,nrow(pdiag))))
			})
	deltamat=tapply(deltas,deltaddl$time,
			function(delta) 
			{
				deltamat=matrix(1,ncol=length(delta),nrow=length(delta))
				diag(deltamat)=delta
				return(cbind(rbind(rep(1,ncol(deltamat),deltamat,1-delta),rep(1,nrow(pmat)))))
			})
	return(lapply(1:length(pmat), function(i) pmat[[i]]*deltamat[[i]]))
}
