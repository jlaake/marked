###################################################################################
# cjs.initial - computes starting values for Phi and p parameters from the
# list of design matrices and the summarized data list including ch matrix and
# first and last vectors. If any values are missing (NA) or abs(par)>5, they are
# set to 0.
#
# Arguments:
#    dml             - ist of design matrices for Phi and p
#    imat            - list containing chmat, first and last 
#
# Value: initial parameter vector
#
####################################################################################
cjs.initial=function(dml,imat)
{
#   Create initial values for p using bernoulli glm with Manley-Parr approach
	cat("\n Computing initial parameter estimates\n")
	num=nrow(imat$chmat)
	ind=matrix(c(1:num,imat$first+1,imat$last-1),ncol=3,nrow=num)
	ind=ind[ind[,2]<=ind[,3],]
	dep.values=unlist(apply(ind,1,function(x,z)z[x[1],x[2]:x[3]],z=imat$chmat))
	dd.indices=unlist(apply(ind,1,function(x) (x[1]-1)*(ncol(imat$chmat)-1)+x[2]:x[3]))-1
	x=dml[["p"]][dd.indices,]
    initial.p=coef(glm.fit(x,dep.values,family=binomial()))
	initial.p[is.na(initial.p)]=0
	initial.p[initial.p < -5 | initial.p > 5]=0
	#   Create initial values for S using bernoulli glm assuming p=1
	last=imat$last+1
	last[last>ncol(imat$chmat)]=ncol(imat$chmat)
	ind=matrix(c(1:num,imat$first+1,last),ncol=3,nrow=num)
	ind=ind[ind[,2]<=ncol(imat$chmat),]
	dep.values=unlist(apply(ind,1,function(x,z)c(rep(1,(x[3]-x[2])),z[x[1],x[3]]),z=imat$chmat))
	dd.indices=unlist(apply(ind,1,function(x) (x[1]-1)*(ncol(imat$chmat)-1)+x[2]:x[3]))-1
	x=dml[["Phi"]][dd.indices,]
	initial.Phi=coef(glm.fit(x,dep.values,family=binomial()))
	initial.Phi[is.na(initial.Phi)]=0
	initial.Phi[initial.Phi < -5 | initial.Phi > 5]=0
	return(c(initial.Phi,initial.p))
}
