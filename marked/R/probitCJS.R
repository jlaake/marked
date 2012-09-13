#' Perform MCMC analysis of a CJS model
#' 
#' Takes design data list created with the function \link{make.design.data} for model "probitCJS" 
#' and draws a sample from the posterior distribution using a Gibbs sampler.
#' 
#' 
#' @param data A list of design data from make.design.data for model="probitCJS"
#' @param parameters A model specification list with a list for Phi and p containing a formula and optionally a prior specification which is a named list
#'        containing 'mu', the prior mean and 'tau' the scale for the conjugate prior precision matrix X'X. 
#' @param burnin number of iteration to initially discard for MCMC burnin
#' @param iter number of iteration to run the Gibbs sampler for following burnin
#' @param init.list A named list (beta.z, beta.y). If null and imat is not null, uses cjs.initial to create initial values; otherwise assigns 0
#' @param imat A list of vectors and matrices constructed by \code{\link{process.ch}} from the capture history data
#' @return A list with MCMC iterations and summarized output:
#' \item{beta.mcmc}{list with elements Phi and p which contain MCMC iterations for each beta parameter} 
#' \item{real.mcmc}{list with elements Phi and p which contain MCMC iterations for each real parameter} 
#' \item{beta}{list with elements Phi and p which contain summary of MCMC iterations for each beta parameter} 
#' \item{reals}{list with elements Phi and p which contain summary of MCMC iterations for each real parameter} 
#' @export
#' @author Devin Johnson
#' @examples
#'
#' # Analysis of the dipper data
#' data(dipper)
#' # following example uses unrealistically low values for burnin and iteration to reduce package testing time
#' fit1 <- crm(dipper,model="probitCJS",model.parameters=list(Phi=list(formula=~time*sex),p=list(formula=~time+sex)), burnin=50, iter=250)
#' fit1
#' # Real parameter summary
#' fit1$reals
probitCJS <- function(data, parameters=list(Phi=list(formula=~1),p=list(formula=~1)), burnin, iter, init.list=NULL, imat=NULL){
  
  ### DEFINE SOME FUNCTIONS ###
  
  sample.z <- function(id, mu.y, mu.z, yvec){
  	.Call("sampleZ", ID=id, PVec=pnorm(mu.y), PhiVec=pnorm(mu.z), yVec=yvec, PACKAGE="marked")
  	}
  make.ztilde.idx <- function(id, zvec){
  	.Call("makeZtildeIdx", ID=id, zvec=zvec, PACKAGE="marked")
  	}

  #### changed 12 July 2012 to use cjs.initial if init.list is null and imat available
  ### INITIAL VALUES ### changed to NULL from missing
  phi.model <- parameters$Phi$formula
  p.model <- parameters$p$formula
  if(is.null(init.list)){
	  Xy <- model.matrix(p.model, data$p)
	  Xz <- model.matrix(phi.model, data$Phi)	  
	  if(is.null(imat)){
		  beta.z <- rep(0,ncol(Xz))
		  beta.y <- rep(0,ncol(Xy))	
	  }else
	  {
		  beta <- cjs.initial(list(Phi=Xz,p=Xy),imat=imat,link="probit")
		  beta.z <- beta[1:ncol(Xz)]	
		  beta.y <- beta[(ncol(Xz)+1):length(beta)]
	  }
  }
  else{
	  beta.z <- init.list$beta.z
	  beta.y <- init.list$beta.y
  }
  
  ### DATA MANIPULATION ###
  #### added by jll 2 July 2012 so it uses design data (ie make.design.data)
	p.data <- data$p[data$p$Time>=data$p$Cohort,]
  p.data <- p.data[order(p.data$id, p.data$Time),]
	phi.data <- data$Phi[data$Phi$Time>=data$Phi$Cohort,]
  phi.data <- phi.data[order(phi.data$id, phi.data$Time),]
	yvec <- p.data$Y
	n <- length(yvec)
	Xy <- model.matrix(p.model, p.data)
	pn.p <- colnames(Xy)
	Xz <- model.matrix(phi.model, phi.data)
	pn.phi <- colnames(Xz)
	id <- p.data$id
	uXy <- unique(Xy)
	colnames(uXy) <- pn.p
	uXz <- unique(Xz)
	colnames(uXz) <- pn.phi
  
  ###  PRIOR DISTRIBUTIONS ### - changed jll 2 July 2012
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
	  
  ### STORAGE ###
  beta.z.stor <- matrix(NA, iter, ncol(Xz))
  beta.y.stor <- matrix(NA, iter, ncol(Xy))
  fitted.y.stor <- matrix(NA, iter, nrow(uXy))
  fitted.z.stor <- matrix(NA, iter, nrow(uXz))
  
  colnames(beta.z.stor) <- pn.phi
  colnames(beta.y.stor) <- pn.p
  
  ### BEGIN MCMC ###
  cat("probitCJS MCMC beginning...\n")
  cat("p model = ", as.character(p.model),"\n")
  cat("phi model = ", as.character(phi.model),"\n")
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
    if(m>burnin){
      beta.z.stor[m-burnin,] <- beta.z
      fitted.z.stor[m-burnin,] <- pnorm(uXz%*%beta.z)
    }
    
    ### UPDATE Y.TILDE ### 
    a <- ifelse(yvec==0, -Inf, 0)
    b <- ifelse(yvec==0, 0, Inf)
    y.tilde <- rtruncnorm(n, a=a, b=b, mean=Xy%*%beta.y, sd=1)
    
    ### BETA.Y UPDATE ###
    Q.b.y <- tau.b.y*crossprod(Xy[zvec==1,])
    V.beta.y.inv <- crossprod(Xy[zvec==1,]) + Q.b.y
    m.beta.y <- solve(V.beta.y.inv, crossprod(Xy[zvec==1,],y.tilde[zvec==1])+crossprod(Q.b.y,mu.b.y))
    beta.y <- m.beta.y + solve(chol(V.beta.y.inv), rnorm(ncol(Xy),0,1))
    if(m>burnin) {
      beta.y.stor[m-burnin,] <- beta.y
      fitted.y.stor[m-burnin,] <- pnorm(uXy%*%beta.y)
    }
    
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
  
  mode <- function(x){
  	dx <- density(x)
  	dx$x[dx$y==max(dx$y)]
  	}
  
  fitted.phi <- mcmc(fitted.z.stor)
  summ.phi <- summary(fitted.phi)
  hpd.phi <- HPDinterval(fitted.phi)
  uphimf <- unique(model.frame(phi.model, phi.data))
  if(phi.model==as.formula(~1)) uphimf <- data.frame(Intercept=1)
  reals.phi <- data.frame(uphimf, 
		  mode = apply(fitted.z.stor, 2, mode), 
		  mean=apply(fitted.z.stor, 2, mean), 
      sd = apply(fitted.z.stor, 2, sd),
		  CI.lower=hpd.phi[,1], CI.upper=hpd.phi[,2])
  reals.phi <- reals.phi[do.call(order, reals.phi),]
  rownames(reals.phi) <- NULL
  fitted.p <- mcmc(fitted.y.stor)
  summ.p <- summary(fitted.p)
  hpd.p <- HPDinterval(fitted.p)
  upmf <- unique(model.frame(p.model, p.data))
  if(p.model==as.formula(~1)) upmf <- data.frame(Intercept=1)
  reals.p <- data.frame(upmf, 
		  mode = apply(fitted.y.stor, 2, mode), 
		  mean=apply(fitted.y.stor,2,mean), 
      sd=apply(fitted.y.stor,2,sd),
		  CI.lower=hpd.p[,1], CI.upper=hpd.p[,2])
  reals.p <- reals.p[do.call(order, reals.p),]
  rownames(reals.p) <- NULL
  phibeta.mcmc <- mcmc(beta.z.stor)
  summ.phi <- summary(phibeta.mcmc)
  hpd.phi <- HPDinterval(phibeta.mcmc)
  beta.phi <- data.frame(posterior.mean=apply(beta.z.stor, 2, mean), 
		  CI.lower=hpd.phi[,1], CI.upper=hpd.phi[,2])
  
  pbeta.mcmc <- mcmc(beta.y.stor)
  summ.p <- summary(pbeta.mcmc)
  hpd.p <- HPDinterval(pbeta.mcmc)
  beta.p <- data.frame(  posterior.mean=apply(beta.y.stor, 2, mean), 
		  CI.lower=hpd.p[,1], CI.upper=hpd.p[,2])
  
  res=list(beta.mcmc=list(Phi= phibeta.mcmc,p= pbeta.mcmc), 
		   reals.mcmc=list(Phi=fitted.phi,p=fitted.p),
           beta=list(Phi=beta.phi,p=beta.p),
		   reals=list(Phi=reals.phi,p=reals.p))
  class(res)=c("crm","probitCJS")
  return(res)
}	### END OF FUNCTION ###
