cjs.lnl=function(par,imat,Phi.dm,p.dm,Phi.fixed=NULL,p.fixed=NULL,Phi.links=NULL,p.links=NULL,debug=FALSE,
		time.intervals=NULL,all=FALSE)
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
	nocc=imat$nocc
	nphi=ncol(Phi.dm)
	np=ncol(p.dm)
	beta.phi=par[1:nphi]
	beta.p=par[(nphi+1):(nphi+np)]
	Phibeta=as.vector(Phi.dm%*%beta.phi)
	if(length(Phi.links>0))
	{
		Phibeta[Phi.links]=as.vector((sin(Phi.dm[Phi.links,,drop=FALSE]%*%beta.phi)+1)/2)
	    Phibeta[Phi.links]=log(1/(1/Phibeta[Phi.links]-1))
	}
	Phibeta=matrix(Phibeta,ncol=nocc-1,nrow=nrow(Phi.dm)/(nocc-1),byrow=TRUE)
	pbeta=as.vector(p.dm%*%beta.p)
	if(length(p.links>0))
	{
	   pbeta[p.links]=as.vector((sin(p.dm[p.links,,drop=FALSE]%*%beta.p)+1)/2)
	   pbeta[p.links]=log(1/(1/pbeta[p.links]-1))
    }
	pbeta=matrix(pbeta,ncol=nocc-1,nrow=nrow(p.dm)/(nocc-1),byrow=TRUE)	
#	pbeta=matrix(as.vector(p.dm%*%beta.p),ncol=nocc-1,nrow=nrow(p.dm)/(nocc-1),byrow=TRUE)
	if(is.null(time.intervals))time.intervals=matrix(1,nrow=length(imat$first),ncol=nocc-1)
	lnl=.Fortran("cjs",as.double(imat$chmat),as.double(Phibeta),as.double(pbeta),
			as.double(imat$first),as.double(imat$last),as.double(imat$freq),
			as.integer(imat$loc),as.double(Phi.fixed),as.double(p.fixed),
			as.double(time.intervals),as.integer(nrow(imat$chmat)),           
			as.integer(ncol(imat$chmat)),as.integer(nrow(Phi.fixed)),
			as.integer(nrow(p.fixed)),lnl=double(1),p0=double(nrow(imat$chmat)),PACKAGE="marked")
	if(debug)
	{
		cat("\nlnl = ",lnl$lnl)
	} else
	if((f_eval-100*floor(f_eval/100))==0)
	{
		string=paste(" Number of evaluations: ",f_eval," -lnl:",formatC(lnl$lnl,digits=10))
		cat(paste(paste(rep("\n",nchar(string)),collapse=""),string,sep=""))
		flush.console()
	}	
	if(all)
		return(lnl)
	else
		return(lnl$lnl)
}
