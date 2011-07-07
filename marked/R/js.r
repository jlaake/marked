js=function(x,ddl,dml,parameters,accumulate=TRUE,Phi=NULL,p=NULL,initial=NULL,method="BFGS",
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,...)
###################################################################################
# js - convenience function for optim to optimize the js likelihood (js.lnl) for
#       a given model specified by the Phi.dm, p.dm, pent.dm, and N.dm design matrices and
#       the data x.
#
# Arguments:
#    x               - processed dataframe created by process.data; data contains 
#                         field ch and freq at minimum
#    ddl             - list of dataframes for design data
#    dml             - list of design matrices
#    parameters      - list of parameter model specifications
#    accumulate     - if TRUE will accumulate capture histories with common value
#                      and with a common design matrix and fixed values for Phi and p
#    Phi             - initial value for intercept 
#    p               - initial value for intercept
#    initial         - initial values for parameters if desired; if it is a
#                      named vector it will match with column names of design matrices
#    method          - method used in optim
#    hessian         - if TRUE, computes and returns hessian
#    debug           - if TRUE show iterations with par and -2lnl
#    chunk_size      - specifies amount of memory to use in accumulating capture histories
#                             use is 8*chunk_size/1e6 MB (default 80MB)
#  refit             - number of times to refit the model if it doesn't converge
#    ...             - additional parameters passed to js.lnl or optim
#
# Value: list of results
#
###################################################################################
{
#
#  Setup values from arguments
#    Phi.dm          - Phi design matrix created by call to create.dm
#    p.dm            - p design matrix created by call to create.dm
#    pent.dm         - pent design matrix created by call to create.dm
#    N.dm            - N design matrix created by call to create.dm
#    Phi.dmdf        - Phi design dataframe created by call to create.dmdf
#    p.dmdf          - p design dataframe created by call to create.dmdf
#    pent.dmdf       - pent design dataframe created by call to create.dmdf
#    N.dmdf          - N design dataframe created by call to create.dmdf
#    Phi.fixed       - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix phi(i,j)=f 
#    p.fixed         - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
#    pent.fixed      - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix pent(i,j)=f 
 Phi.dmdf=ddl$Phi
 p.dmdf=ddl$p
 Phi.dm=dml$Phi
 p.dm=dml$p
 pent.dmdf=ddl$pent
 N.dmdf=ddl$N
 pent.dm=dml$pent
 N.dm=dml$N
 nocc=x$nocc
 #  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal
 time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
 if(!is.null(Phi.dmdf$time.interval))
	 time.intervals=matrix(Phi.dmdf$time.interval,nrow(x$data),ncol=nocc-1,byrow=T)
 Phi.fixed=parameters$Phi$fixed
 p.fixed=parameters$p$fixed
 pent.fixed=parameters$pent$fixed
 N.fixed=parameters$N$fixed
 if(nrow(N.dm)==1) 
    nobstot=sum(x$data$freq)    
 else
   nobstot=tapply(x$data$freq,x$data$group,sum)   
   x=x$data
   if(is.null(Phi.fixed))
      Phi.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)
   if(is.null(p.fixed))
      p.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
   if(is.null(pent.fixed))
      pent.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
   freq=x$freq
#  if data are to be accumulated based on ch and design matrices do so here;
#  otherwise simply set ch
   if(accumulate)
   {
	  cat("\n Accumulating capture frequencies based on design. This can take awhile.\n")
	  flush.console()
#     Save dms, fixed values and imat before accumulating
      Phi.dm.save=Phi.dm
      p.dm.save=p.dm
      pent.dm.save=pent.dm
      N.dm.save=N.dm
	  Phi.fixed.save=Phi.fixed
	  p.fixed.save=p.fixed
	  pent.fixed.save=pent.fixed
#     Create chmat and ch
	  ch=x$ch
	  imat=process.ch(ch,freq,all=TRUE)
	  imat.save=imat
      chmat=imat$chmat
      if(sum(imat$loc)>0)
      {
         indices=which(imat$loc==1)
         chmat[cbind(indices,imat$last[indices])]=2
         ch=apply(chmat,1,paste,collapse="")
      }
#      else
#         ch=x$ch
#     Create fixed parameter matrices
	  if(Phi.fixed[1,1]>0)
      {
         Phifixmat=rep(NA,nrow(x)*(nocc-1))
         Phifixmat[(nocc-1)*(Phi.fixed[,1]-1)+Phi.fixed[,2]]=Phi.fixed[,3]      
      }
      else
		  Phifixmat=NULL
	  if(p.fixed[1,1]>0)
      {         
         pfixmat=rep(NA,nrow(x)*nocc)
         pfixmat[nocc*(p.fixed[,1]-1)+p.fixed[,2]-1]=p.fixed[,3]      
      }else
		 pfixmat=NULL
      if(pent.fixed[1,1]>0)
      {         
         pentfixmat=rep(NA,nrow(x)*(nocc-1))
         pentfixmat[(nocc-1)*(pent.fixed[,1]-1)+pent.fixed[,2]-1]=pent.fixed[,3]      
      }else
		  pentfixmat=NULL
#     Start accumulating 
#     Construct common splits based on Phi and pent dms and fixed parameter values
	  nrows=nrow(Phi.dm)
	  pieces=floor(ncol(Phi.dm)*nrows/chunk_size)+1
	  chdesign=NULL
	  iseq=seq(1,nrow(x),floor(nrows/pieces/(nocc-1)))
	  if(iseq[length(iseq)]!=nrow(x))iseq=c(iseq,nrow(x))
	  prev_i=1
	  for(i in iseq[-1])
	  {
		  lower=(prev_i-1)*(nocc-1)+1
		  upper=i*(nocc-1)
		  xdesign=cbind(rep(ch[prev_i:i],each=nocc-1),as.matrix(Phi.dm[lower:upper,,drop=FALSE]),as.matrix(pent.dm[lower:upper,,drop=FALSE]),
				  Phifixmat[lower:upper],pentfixmat[lower:upper])
		  xdesign=sapply(split(xdesign,rep(1:(i-prev_i+1),each=nocc-1)),paste,collapse="")
		  chdesign=c(chdesign,xdesign)
		  prev_i=i+1
	  }
#     If time intervals vary across individuals, then use it as part of accumulation 
#     split
	  occ.int=as.vector(unlist(apply(time.intervals,2,unique)))
	  if(length(occ.int)!=(nocc-1))
		  chdesign=paste(chdesign,apply(time.intervals,1,paste,collapse=""),sep="")
#     Now use p which has different parameter structure	 	  
      nrows=nrow(p.dm)
      pieces=floor(ncol(p.dm)*nrows/chunk_size)+1
      chdesignp=NULL
      iseq=seq(1,nrow(x),floor(nrows/pieces/(nocc)))
      if(iseq[length(iseq)]!=nrow(x))iseq=c(iseq,nrow(x))
      prev_i=1
      for(i in iseq[-1])
      {
	     lower=(prev_i-1)*nocc+1
	     upper=i*nocc
	     xdesign=cbind(rep(ch[prev_i:i],each=nocc),as.matrix(p.dm[lower:upper,,drop=FALSE]),
		  	   pfixmat[lower:upper])
	     xdesign=sapply(split(xdesign,rep(1:(i-prev_i+1),each=nocc)),paste,collapse="")
	     chdesignp=c(chdesignp,xdesign)
	     prev_i=i+1
      }
#     Split data based on paste of all dms, fixed values
	  chsplit=split(1:nrow(x),paste(chdesign,chdesignp,as.numeric(freq==0),sep=""))
#     Get new set of indices
	  indices=as.vector(sapply(chsplit,min))
#     Accumulate frequencies
	  counts=as.vector(sapply(chsplit,function(x)sum(freq[x])))
      freq=counts[order(indices)]
#     Get reduced set of dms and fixed parameter matrices
	  dm.index=rep(1:nrow(x),each=nocc-1)%in%indices
      Phi.dm=Phi.dm[dm.index,,drop=FALSE]
      pent.dm=pent.dm[dm.index,,drop=FALSE]
      if(!is.null(Phifixmat))
      {
        chi=(rep(1:nrow(x),each=nocc-1))[dm.index][!is.na(Phifixmat)[dm.index]]
        occj=(rep(1:(nocc-1),times=nrow(x)))[dm.index][!is.na(Phifixmat)[dm.index]]
        val=Phifixmat[dm.index][!is.na(Phifixmat)[dm.index]]
        Phi.fixed=cbind( match(chi,indices),occj,val)
       }  
      dm.index=rep(1:nrow(x),each=nocc)%in%indices
      p.dm=p.dm[dm.index,,drop=FALSE]
      if(!is.null(pfixmat))
      {
        chi=(rep(1:nrow(x),each=nocc))[dm.index][!is.na(pfixmat)[dm.index]]
        occj=(rep(1:nocc,times=nrow(x)))[dm.index][!is.na(pfixmat)[dm.index]]
        val=pfixmat[dm.index][!is.na(pfixmat)[dm.index]]
        p.fixed=cbind( match(chi,indices),occj,val)
       }  
      ch=x$ch[sort(indices)]
#     Compute new imat
	  imat=process.ch(ch,freq,all=TRUE)
	  nx=sum(x$freq)
	  if(sum(freq)!=nx)
		   stop(paste("Error in accumulation. Number of accumulated",sum(freq),"not equal to original number",nx))
	  cat(" ",nrow(x)," capture histories collapsed into ",length(ch),"\n")
	  
   }else
   {
#  get first and last vectors, loc and chmat
	   ch=x$ch
	   imat=process.ch(ch,freq,all=TRUE)
	   imat.save=NULL
	   Phi.dm.save=NULL
	   p.dm.save=NULL  
	   pent.dm.save=NULL  
	   N.dm.save=NULL  
	   Phi.fixed.save=NULL
	   p.fixed.save=NULL  
	   pent.fixed.save=NULL  
   }
#  if no initial values set, set some default ones
   if(is.null(initial))
   {
      if(is.null(Phi))Phi=0.5
      if(is.null(p))p=0.5   
   }           
#  call optim to find mles with js.lnl which gives -2 * log-likelihood
   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
   cat("\n Starting optimization",ncol(Phi.dm)+ncol(p.dm)+ncol(pent.dm)+ncol(N.dm)," parameters\n")
   if(is.null(initial))
      par=c(log(Phi/(1-Phi)),rep(0,ncol(Phi.dm)-1),log(p/(1-p)),rep(0,ncol(p.dm)-1),
        0,rep(0,ncol(pent.dm)-1),1,rep(0,ncol(N.dm)-1))
   else
   {
      if(is.null(names(initial)))
      {
         if(length(initial)!=(ncol(Phi.dm)+ncol(p.dm)+ncol(pent.dm)+ncol(N.dm)))
            stop("Length of initial vector does not match number of parameters.")
         else
            par=initial
      }
      else
      {
        beta.names=c(paste("Phi:",colnames(Phi.dm),sep="") ,paste("p:",colnames(p.dm),sep=""))
        par=rep(0,length(beta.names))
        par[beta.names%in%names(initial)]=initial[which(names(initial)%in%beta.names)]
      }
   }
   convergence=1
   i=0
   while (convergence!=0 & i <= refit)
   {
	   mod=optimx(par,js.lnl,imat=imat,Phi.dm=Phi.dm,p.dm=p.dm,pent.dm=pent.dm,N.dm=N.dm,Phi.fixed=Phi.fixed,
			   p.fixed=p.fixed,pent.fixed=pent.fixed,method=method,hessian=hessian,
			   debug=debug,time.intervals=time.intervals,control=control,itnmax=itnmax,nobstot=nobstot,...)
	   par=mod$par$par
	   convergence=mod$conv
	   i=i+1
   }
   js.beta=mod$par$par
   names(js.beta)=c(paste("Phi:",colnames(Phi.dm),sep="") ,paste("p:",colnames(p.dm),sep=""),paste("pent:",colnames(pent.dm),sep=""),paste("N:",colnames(N.dm),sep=""))
   if(!is.null(imat.save)) imat=imat.save
   if(!is.null(Phi.dm.save)) Phi.dm=Phi.dm.save
   if(!is.null(p.dm.save)) p.dm=p.dm.save
   if(!is.null(pent.dm.save)) pent.dm=pent.dm.save
   if(!is.null(Phi.fixed.save)) Phi.fixed=Phi.fixed.save
   if(!is.null(p.fixed.save)) p.fixed=p.fixed.save
   if(!is.null(pent.fixed.save)) pent.fixed=pent.fixed.save
   nphi=ncol(Phi.dm)
   np=ncol(p.dm)
   npent=ncol(pent.dm)
   nN=ncol(N.dm)
   beta.phi=js.beta[1:nphi]
   beta.p=js.beta[(nphi+1):(nphi+np)]
   beta.pent=js.beta[(nphi+np+1):(nphi+np+npent)]
   beta.N=js.beta[(nphi+np+npent+1):(nphi+np+npent+nN)]
# create Phi and p beta matrices excluding first occasion on p
   Phis=plogis(as.vector(Phi.dm%*%beta.phi))
   ps=plogis(as.vector(p.dm%*%beta.p))
   pents=cbind(rep(1,nrow(pent.dm)),exp(as.vector(pent.dm%*%beta.pent)))
   pents=pents/apply(pents,1,sum)
   Ns=exp(as.vector(N.dm%*%beta.N))
   reals=p.dmdf
   names(reals)[names(reals)=="time"]="p.time"
   names(reals)[names(reals)=="age"]="p.age"
   Phirecs=Phi.dmdf[,!colnames(Phi.dmdf)%in%colnames(reals),drop=FALSE]
   Phirecs[(rep(1:nrow(x),each=nocc)-1)*(nocc-1)+c(1:(nocc-1),nocc-1),]
   Phirecs[seq(nocc,nrow(x)*nocc,by=nocc),]=NA
#   split(Phi.dmdf[,!colnames(Phi.dmdf)%in%colnames(reals),drop=FALSE],rep(1:nrow(x),each=nocc-1))
#   reals=cbind(reals,
   names(reals)[names(reals)=="time"]="Phi.time"
   names(reals)[names(reals)=="age"]="Phi.age"
#   reals$Phi=as.vector(t(matrix(allval[[2]],nrow=nrow(x),ncol=nocc-1)))
#   reals$p=as.vector(t(matrix(allval[[3]],nrow=nrow(x),ncol=nocc-1)))
   ui=tapply(imat$freq,list(imat$first,x$group),sum)
#   ui=tapply(imat$freq,factor(imat$first),sum)
   lnl=mod$fvalues$fvalues+ 2*sum(lfactorial(ui))
   res=list(beta=js.beta,lnl=lnl,AIC=lnl+2*length(js.beta),convergence=mod$conv,count=mod$itns,reals=reals)
   if(!is.null(mod$hessian)) 
   {
      res$vcv=solve(mod$hessian)
      colnames(res$vcv)=names(js.beta)
      rownames(res$vcv)=names(js.beta)
   }   
   class(res)=c("crm","js")
   return(res)
}


