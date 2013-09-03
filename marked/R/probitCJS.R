#' Perform MCMC analysis of a CJS model
#' 
#' Takes design data list created with the function \link{make.design.data} for model "probitCJS" 
#' and draws a sample from the posterior distribution using a Gibbs sampler.
#' 
#' 
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param parameters A model specification list with a list for Phi and p containing a formula and optionally a prior specification which is a named list
#'        containing 'mu', the prior mean and 'tau' the scale for the conjugate prior precision matrix X'X. 
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param burnin number of iteration to initially discard for MCMC burnin
#' @param iter number of iteration to run the Gibbs sampler for following burnin
#' @param initial A named list (Phi,p). If null and imat is not null, uses cjs.initial to create initial values; otherwise assigns 0
#' @param imat A list of vectors and matrices constructed by \code{\link{process.ch}} from the capture history data
#' @return A list with MCMC iterations and summarized output:
#' \item{beta.mcmc}{list with elements Phi and p which contain MCMC iterations for each beta parameter} 
#' \item{real.mcmc}{list with elements Phi and p which contain MCMC iterations for each real parameter} 
#' \item{beta}{list with elements Phi and p which contain summary of MCMC iterations for each beta parameter} 
#' \item{reals}{list with elements Phi and p which contain summary of MCMC iterations for each real parameter} 
#' @export
#' @import coda truncnorm
#' @author Devin Johnson
#' @examples
#' \donttest{
#' # Analysis of the dipper data
#' data(dipper)
#' # following example uses unrealistically low values for burnin and iteration to reduce package testing time
#' fit1 <- crm(dipper,model="probitCJS",model.parameters=list(Phi=list(formula=~time*sex),p=list(formula=~time+sex)), burnin=100, iter=1000)
#' fit1
#' # Real parameter summary
#' fit1$results$reals
#' }
probitCJS <- function(ddl,dml,parameters,design.parameters,burnin, iter, initial=NULL, imat=NULL){
  
  ### DEFINE SOME FUNCTIONS ###
  
  sample.z <- function(id, mu.y, mu.z, yvec){
  	.Call("sampleZ", ID=id, PVec=pnorm(mu.y), PhiVec=pnorm(mu.z), yVec=yvec, PACKAGE="marked")
  	}
  make.ztilde.idx <- function(id, zvec){
  	.Call("makeZtildeIdx", ID=id, zvec=zvec, PACKAGE="marked")
  	}
  ### Initial values
  if(is.null(initial)) beta <- cjs.initial(dml,imat=imat,link="probit")
  else beta <- set.initial(c("Phi","p"),dml,initial)$par
  beta.z <- beta$Phi	
  beta.y <- beta$p	  
  ### DATA MANIPULATION ###
  ##  restrict design data so Time>=Cohort and recreate design matrices 
  restricted.dml <- create.dml(ddl,parameters,design.parameters,restrict=TRUE)  
  ddl$p <- ddl$p[ddl$p$Time>=ddl$p$Cohort,]
  yvec <- ddl$p$Y
  n <- length(yvec)
  Xy <- as.matrix(restricted.dml$p$fe)
  pn.p <- colnames(Xy)
  Xz <- as.matrix(restricted.dml$Phi$fe)
  pn.phi <- colnames(Xz)
  id <- ddl$p$id
  ###  PRIOR DISTRIBUTIONS ###
  if(is.null(parameters$Phi$prior)){
	  tau.b.z <- 0.01
	  mu.b.z <- rep(0,ncol(Xz))
  }else{
	  tau.b.z <- parameters$Phi$prior$tau
	  mu.b.z <- parameters$Phi$prior$mu
  }		
  if(is.null(parameters$p$prior)){
	  tau.b.y <- 0.01
	  mu.b.y <- rep(0,ncol(Xy))	
  }else{
	  tau.b.y <- parameters$p$prior$tau
	  mu.b.y <- parameters$p$prior$mu
  }
  if(!is.null(dml$Phi$re)){
    n.Phi.re=length(dml$Phi$re)
    if(is.null(parameters$Phi$priors$re)){
      a.phi=rep(2,n.Phi.re)
      b.phi=rep(1.0E-4,n.Phi.re)
    }
    else{
      a.phi=parameters$Phi$priors$re$a.phi
      b.phi=parameters$Phi$priors$re$b.phi
    }
  }
  if(!is.null(dml$p$re)){
    n.p.re=length(dml$p$re)
    if(is.null(parameters$p$priors$re)){
      a.p=rep(2,n.p.re)
      b.p=rep(1.0E-4,n.p.re)
    }
    else{
      a.p=parameters$p$priors$re$a.p
      b.p=parameters$p$priors$re$b.p
    }
  }
  
	  
  ### STORAGE ###
  beta.z.stor <- matrix(NA, iter, ncol(Xz))
  beta.y.stor <- matrix(NA, iter, ncol(Xy))
  colnames(beta.z.stor) <- pn.phi
  colnames(beta.y.stor) <- pn.p
  
  ### BEGIN MCMC ###
  cat("probitCJS MCMC beginning...\n")
  cat("p model = ", as.character(parameters$p$formula),"\n")
  cat("phi model = ", as.character(parameters$Phi$formula),"\n")
  flush.console()
  
  tot.iter <- burnin + iter
  st <- Sys.time()
  for(m in 1:tot.iter){
    
    ### UPDATE Z ###
   zvec <- sample.z(id=id, mu.y=Xy%*%beta.y, mu.z=Xz%*%beta.z, yvec)
    
    ### UPDATE Z.TILDE ### 
    a <- ifelse(zvec==0, -Inf, 0)
    b <- ifelse(zvec==0, 0, Inf)
    z.tilde <- rtruncnorm(n, a=a, b=b, mean=Xz%*%beta.z, sd=1)
    
    ### BETA.Z UPDATE ###
    idx.z.tilde <- make.ztilde.idx(id, zvec)
    Q.b.z <- tau.b.z*crossprod(Xz[idx.z.tilde,])
    V.beta.z.inv <- crossprod(Xz[idx.z.tilde,]) + Q.b.z
    m.beta.z <- solve(V.beta.z.inv, crossprod(Xz[idx.z.tilde,],z.tilde[idx.z.tilde]) + crossprod(Q.b.z,mu.b.z))
    beta.z <- m.beta.z + solve(chol(V.beta.z.inv), rnorm(ncol(Xz),0,1))
    if(m>burnin)beta.z.stor[m-burnin,] <- beta.z
    
    ### UPDATE Y.TILDE ### 
    a <- ifelse(yvec==0, -Inf, 0)
    b <- ifelse(yvec==0, 0, Inf)
    y.tilde <- rtruncnorm(n, a=a, b=b, mean=Xy%*%beta.y, sd=1)
    
    ### BETA.Y UPDATE ###
    Q.b.y <- tau.b.y*crossprod(Xy[zvec==1,])
    V.beta.y.inv <- crossprod(Xy[zvec==1,]) + Q.b.y
    m.beta.y <- solve(V.beta.y.inv, crossprod(Xy[zvec==1,],y.tilde[zvec==1])+crossprod(Q.b.y,mu.b.y))
    beta.y <- m.beta.y + solve(chol(V.beta.y.inv), rnorm(ncol(Xy),0,1))
    if(m>burnin) beta.y.stor[m-burnin,] <- beta.y
   
   ### RANDOM EFFECT UPDATES ###
   
    
    ### TIMING OF SAMPLER ###
    if(m==30){
      tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/30
      ttc <- round((tot.iter-30)*tpi/3600, 2)
      if(ttc>=1) cat("Approximate time till completion: ", ttc, " hours\n")
      else cat("Approximate time till completion: ", ttc*60, " minutes\n")
    }
    if(100*(m/tot.iter) >= 10 & (100*(m/tot.iter))%%10==0) 
	{
		cat(100*(m/tot.iter), "% completed\n")
		flush.console()
	}
  }
  ### MAKE SOME OUTPUT ### jll 2 July changed to summarize betas as well 
  ### 21 Sept jll removed real computations; it is now in compute.real called from crm
  phibeta.mcmc <- mcmc(beta.z.stor)
  summ.phi <- summary(phibeta.mcmc)
  hpd.phi <- HPDinterval(phibeta.mcmc)
  beta.phi <- data.frame(mode = apply(beta.z.stor, 2, mcmc_mode), mean=apply(beta.z.stor, 2, mean), 
		  sd=apply(beta.z.stor,2,sd),CI.lower=hpd.phi[,1], CI.upper=hpd.phi[,2])
  pbeta.mcmc <- mcmc(beta.y.stor)
  summ.p <- summary(pbeta.mcmc)
  hpd.p <- HPDinterval(pbeta.mcmc)
  beta.p <- data.frame( mode = apply(beta.y.stor, 2, mcmc_mode), 
		  mean=apply(beta.y.stor, 2, mean), 
		  sd=apply(beta.y.stor,2,sd),
		  CI.lower=hpd.p[,1], CI.upper=hpd.p[,2])  
  res=list(beta.mcmc=list(Phi= phibeta.mcmc,p= pbeta.mcmc), 
		   beta=list(Phi=beta.phi,p=beta.p),
		   model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe))
  class(res)=c("crm","mcmc","probitCJS")
  return(res)
}	### END OF FUNCTION ###
