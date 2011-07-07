cjs=function(x,ddl,dml,parameters,accumulate=TRUE,Phi=NULL,p=NULL,initial=NULL,method,
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale, ...)
###################################################################################
# cjs - convenience function for optim to optimize the cjs likelihood (cjs.lnl) for
#       a given model specified by the Phi.dm and p.dm design matrices and
#       the data x.
#
# Arguments:
#    x               - processed dataframe created by process.data
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
#    refit           - number of times to refit the model if it fails to converge
#    ...             - additional arguments passed to cjs.lnl or optim
#
# Value: list of results
#
###################################################################################
{
#
#  Setup values from arguments
#    Phi.dm          - Phi design matrix created by call to create.dm
#    p.dm            - p design matrix created by call to create.dm
#    Phi.dmdf        - Phi design dataframe created by call to create.dmdf
#    p.dmdf          - p design dataframe created by call to create.dmdf
#    Phi.fixed       - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix phi(i,j)=f 
#    p.fixed         - matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
#                       to fix p(i,j)=f 
   Phi.dmdf=ddl$Phi
   p.dmdf=ddl$p
   Phi.dm=dml$Phi
   p.dm=dml$p
   nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal
   time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
   if(!is.null(Phi.dmdf$time.interval))
	   time.intervals=matrix(Phi.dmdf$time.interval,nrow(x$data),ncol=nocc-1,byrow=T)
   Phi.fixed=parameters$Phi$fixed
   p.fixed=parameters$p$fixed
   if(is.null(Phi.fixed))
      Phi.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)
   if(is.null(p.fixed))
      p.fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)  
   x=x$data
#  set default frequencies if not used
   if(is.null(x$freq))
      freq=NULL
   else
      freq=x$freq
#  get first and last vectors, loc and chmat
   ch=x$ch
   imat=process.ch(ch,freq)
   imat.save=NULL
   Phi.dm.save=NULL
   p.dm.save=NULL  
   Phi.fixed.save=NULL
   p.fixed.save=NULL  
#  if data are to be accumulated based on ch and design matrices do so here;
#  otherwise simply set ch
   if(accumulate)
   {
	  cat("\n Accumulating capture frequencies based on design. This can take awhile.\n")
	  flush.console()
	  imat.save=imat   
      Phi.dm.save=Phi.dm
      p.dm.save=p.dm
      chmat=imat$chmat
      if(sum(imat$loc)>0)
      {
         indices=which(imat$loc==1)
         chmat[cbind(indices,imat$last[indices])]=2
         ch=apply(chmat,1,paste,collapse="")
      }
      else
         ch=x$ch
      Phifixmat=NULL
      pfixmat=NULL
      Phi.fixed.save=Phi.fixed
      p.fixed.save=p.fixed
      if(Phi.fixed[1,1]>0)
      {
         Phifixmat=rep(NA,nrow(x)*(nocc-1))
         Phifixmat[(nocc-1)*(Phi.fixed[,1]-1)+Phi.fixed[,2]]=Phi.fixed[,3]      
      }
      if(p.fixed[1,1]>0)
      {         
         pfixmat=rep(NA,nrow(x)*(nocc-1))
         pfixmat[(nocc-1)*(p.fixed[,1]-1)+p.fixed[,2]-1]=p.fixed[,3]      
      }
	  nrows=nrow(Phi.dm)
	  pieces=floor(max(ncol(Phi.dm),ncol(p.dm))*nrows/chunk_size)+1
	  chdesign=NULL
	  iseq=seq(1,nrow(x),floor(nrows/pieces/(nocc-1)))
      if(iseq[length(iseq)]!=nrow(x))iseq=c(iseq,nrow(x))
	  prev_i=1
	  for(i in iseq[-1])
	  {
 	   lower=(prev_i-1)*(nocc-1)+1
 	   upper=i*(nocc-1)
   	   xdesign=cbind(rep(ch[prev_i:i],each=nocc-1),as.matrix(Phi.dm[lower:upper,,drop=FALSE]),as.matrix(p.dm[lower:upper,,drop=FALSE]),
			                                       Phifixmat[lower:upper],pfixmat[lower:upper])
	   xdesign=sapply(split(xdesign,rep(1:(i-prev_i+1),each=nocc-1)),paste,collapse="")
	   chdesign=c(chdesign,xdesign)
	   prev_i=i+1
	  }
#     If time intervals vary across individuals, then use it as part of accumulation 
#     split
      occ.int=as.vector(unlist(apply(time.intervals,2,unique)))
	  if(length(occ.int)!=(nocc-1))
		  chdesign=paste(chdesign,apply(time.intervals,1,paste,collapse=""),sep="")
      chsplit=split(1:nrow(x),chdesign)
      indices=as.vector(sapply(chsplit,min))
      if(is.null(freq))
         counts=as.vector(sapply(chsplit,length))
      else
         counts=as.vector(sapply(chsplit,function(x)sum(freq[x])))
      freq=counts[order(indices)]
      dm.index=rep(1:nrow(x),each=nocc-1)%in%indices
      Phi.dm=Phi.dm[dm.index,,drop=FALSE]
      p.dm=p.dm[dm.index,,drop=FALSE]
      if(!is.null(Phifixmat))
      {
        chi=(rep(1:nrow(x),each=nocc-1))[dm.index][!is.na(Phifixmat)[dm.index]]
        occj=(rep(1:(nocc-1),times=nrow(x)))[dm.index][!is.na(Phifixmat)[dm.index]]
        val=Phifixmat[dm.index][!is.na(Phifixmat)[dm.index]]
        Phi.fixed=cbind( match(chi,indices),occj,val)
       }  
      if(!is.null(pfixmat))
      {
        chi=(rep(1:nrow(x),each=nocc-1))[dm.index][!is.na(pfixmat)[dm.index]]
        occj=(rep(2:nocc,times=nrow(x)))[dm.index][!is.na(pfixmat)[dm.index]]
        val=pfixmat[dm.index][!is.na(pfixmat)[dm.index]]
        p.fixed=cbind( match(chi,indices),occj,val)
       }  
      ch=x$ch[sort(indices)]
      imat=process.ch(ch,freq) 
	  if(sum(freq)!=nrow(x))stop(paste("Error in accumulation. Number of accumulated",sum(freq),"not equal to original number",nrow(x)))
      cat(" ",nrow(x)," capture histories collapsed into ",length(ch),"\n")
   }
#  Create links   
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
    Phi.links=NULL
	p.links=NULL
#  if no initial values set, set some default ones
   if(is.null(initial))
   {
      if(is.null(Phi))Phi=0.5
      if(is.null(p))p=0.5   
   }     
#  Scale the design matrices with either input scale or computed scale
   if(is.null(scale))
   {
      scale.phi=apply(Phi.dm,2,function(x) mean(x[x!=0]))
	  scale.p=apply(p.dm,2,function(x) mean(x[x!=0]))
   } else
   {
	   scale.phi=scale[1:ncol(Phi.dm)]
	   scale.p=scale[(ncol(Phi.dm)+1):length(scale)]
   }
   for(i in which(scale.phi<.99 | scale.phi>1.01)) Phi.dm[,i]=Phi.dm[,i]/scale.phi[i]
   for(i in which(scale.p<.99 | scale.p>1.01)) p.dm[,i]=p.dm[,i]/scale.p[i]
#  call optim to find mles with cjs.lnl which gives -2 * log-likelihood
   cat("\n Starting optimization for ",ncol(Phi.dm)+ncol(p.dm)," parameters\n")
   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
   flush.console()
   if(is.null(initial))
	  par=c(log(Phi/(1-Phi)),rep(0,ncol(Phi.dm)-1),log(p/(1-p)),rep(0,ncol(p.dm)-1))
   else
   {
      if(is.null(names(initial)))
      {
         if(length(initial)!=(ncol(Phi.dm)+ncol(p.dm)))
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
   if(is.null(scale))par=par*c(scale.phi,scale.p)
   while (convergence!=0 & i <= refit)
   {
      mod=optimx(par,cjs.lnl,imat=imat,Phi.dm=Phi.dm,p.dm=p.dm,Phi.fixed=Phi.fixed,p.fixed=p.fixed,
		   Phi.links=Phi.links,p.links=p.links,method=method,hessian=FALSE,debug=debug,time.intervals=time.intervals,control=control,
		   itnmax=itnmax,...)
	  par=mod$par$par
	  convergence=mod$conv
	  i=i+1
   }
   cjs.beta=mod$par$par/c(scale.phi,scale.p)
   names(cjs.beta)=c(paste("Phi:",colnames(Phi.dm),sep="") ,paste("p:",colnames(p.dm),sep=""))
   if(!is.null(imat.save)) imat=imat.save
   if(!is.null(Phi.dm.save)) Phi.dm=Phi.dm.save
   if(!is.null(p.dm.save)) p.dm=p.dm.save
   if(!is.null(Phi.fixed.save)) Phi.fixed=Phi.fixed.save
   if(!is.null(p.fixed.save)) p.fixed=p.fixed.save
   allval=cjs.lnl(cjs.beta,imat,Phi.dm,p.dm,Phi.fixed,p.fixed,Phi.links,p.links,time.intervals=time.intervals,all=TRUE)
   reals=Phi.dmdf
   names(reals)[names(reals)=="time"]="Phi.time"
   names(reals)[names(reals)=="age"]="Phi.age"
   reals=cbind(reals,p.dmdf[,!colnames(p.dmdf)%in%colnames(reals)])
   names(reals)[names(reals)=="time"]="p.time"
   names(reals)[names(reals)=="age"]="p.age"
   xx=matrix(allval[[2]],nrow=nrow(x),ncol=nocc-1)
   reals$Phi=as.vector(t(xx))
   xx=matrix(allval[[3]],nrow=nrow(x),ncol=nocc-1)
   reals$p=as.vector(t(xx))
   reals=reals[reals$Time>=reals$Cohort,]
   lnl=mod$fvalues$fvalues
   res=list(beta=cjs.beta,lnl=lnl,AIC=lnl+2*length(cjs.beta),convergence=mod$conv,count=mod$itns,reals=reals,mod=mod,scale=c(scale.phi,scale.p))
   if(hessian) 
   {
	   assign(".markedfunc_eval", 0, envir = .GlobalEnv)
	   cat("\n Computing hessian\n")
	   res$hessian=hessian(cjs.lnl,cjs.beta,imat=imat,Phi.dm=Phi.dm,p.dm=p.dm,Phi.fixed=Phi.fixed,p.fixed=p.fixed,Phi.links=Phi.links,
			   p.links=p.links,time.intervals=time.intervals,all=FALSE)
	   res$vcv=solve(res$hessian*outer(1/res$scale,1/res$scale,"*"))
       colnames(res$vcv)=names(cjs.beta)
       rownames(res$vcv)=names(cjs.beta)
   }   
   class(res)=c("crm","cjs")
   return(res)
}


