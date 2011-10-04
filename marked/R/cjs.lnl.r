cjs.lnl=function(par,model_data,Phi.links=NULL,p.links=NULL,debug=FALSE,all=FALSE)
###################################################################################
# cjs.lnl - computes -2*log likelihood value of cjs model using call to FORTRAN
#            subroutine
#
# Arguments:
#    par             - vector of parameter values to evaluate likelihood
#    imat            - list of indicator vectors, frequency etc from process.ch
#    Phi.dm          - Phi design matrix created by call to create.dm
#    p.dm            - p design matrix created by call to create.dm
#    Phi.fixed       - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix phi(i,j)=f 
#    p.fixed         - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
#    debug           - if TRUE show iterations with par and -2lnl
#    time.intervals  - intervals of time between occasions
#    all             - if TRUE, returns entire lnl computation list
#
# Value: -2*log likelihood
##################################################################################
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
		string=paste(" Number of evaluations: ",f_eval," -2lnl:",formatC(lnl$lnl,digits=10))
		cat(paste(paste(rep("\n",nchar(string)),collapse=""),string,sep=""))
		flush.console()
	}	
	if(all)
		return(lnl)
	else
		return(lnl$lnl)
}
