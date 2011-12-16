#' Perform MCMC analysis of a CJS model
#' 
#' Takes an object of class probitCJS.data created with the function \link{make.probitCJS.data} and draws a sample from the posterior distribution using a Gibbs sampler
#' 
#' 
#' @param p.model A formula for the detection portion of the CJS model
#' @param phi.model A formula for the survival portion of the CJS model
#' @param data A probitCJS.data object containg the data for the model
#' @param prior.list A named list containing 
#' \item{phi}{A named list containg 'mu', the prior mean and 'Q' the prior precison matrix for the survival parameters}
#' \item{p}{A named list containg 'mu', the prior mean and 'Q' the prior precison matrix for the survival parameters}
#' @param burnin number of iteration to initially discard for MCMC burnin
#' @param iter number of iteration to run the Gibbs sampler for following burnin
#' @param init.list A named list with initial values [more on this later]
#' @return A list with MCMC output [I'll add to this later] 
#' @export
#' @author Devin Johnson
#' @examples
#'
#' # Analysis of the dipper data
#' data(dipper)
#' dipper <- splitCH(data=dipper)
#' dipper.mcmc <- make.probitCJS.data(dipper, ch.cols=3:9, covariate.cols=2)
#' fit1 <- probitCJS(p.model=~time*sex, phi.model=~time*sex, data=dipper.mcmc, burnin=100, iter=1000)
#' # Real parameter values
#' fit1$reals.phi
#' fit1$reals.p
probitCJS <- function(p.model = ~1, phi.model = ~1, data, prior.list, burnin, iter, init.list){
  
  ### CHECK THE DATA ###
  if(!inherits(data,"probitCJS.data")) stop("The data is not of class 'probitCJS.data'\n
                                           Please see the 'make.probitCJS.data' function.")
### LOAD SOME PACKAGES ###
  #require(Rcpp)
  #require(inline)
  require(truncnorm)
  require(coda)
  
  ### DEFINE SOME FUNCTIONS ###
  
  ln.Phi.1 <- function(mu){pnorm(0, mu, 1, lower.tail=FALSE, log.p=TRUE)}
  ln.Phi.0 <- function(mu){pnorm(0, mu, 1, log.p=TRUE)}
  sample.z <- function(idx, mu.y, mu.z){
    q.y <- ln.Phi.0(mu.y[idx])
    q1.z <- ln.Phi.1(mu.z[idx])
    q0.z <- ln.Phi.0(mu.z[idx])
    pr.z <- exp(c(0,cumsum(q.y)) + c(0,cumsum(q1.z)) + c(q0.z,0))
    d <- sample.int(length(pr.z),1, prob=pr.z)-1
    c(rep(1,d), rep(0, length(pr.z)-1-d))
  }
  
  ### DATA MANIPULATION ###
  
  yvec <- data$Y
  n <- length(yvec)
  zvec.master <- data$Z
  zvec <- zvec.master
  idx.z <- is.na(zvec.master)
  n.unk.z <- sum(idx.z)
  Xy <- model.matrix(p.model, data)
  pn.p <- colnames(Xy)
  Xz <- model.matrix(phi.model, data)
  pn.phi <- colnames(Xz)
  id <- factor(data[,attr(data,"id")])
  id.num <- as.numeric(id)
  id.z <- id[idx.z]
  Xz.z <- matrix(Xz[idx.z,],ncol=ncol(Xz))
  Xy.z <- matrix(Xy[idx.z,],ncol=ncol(Xy))
  uXy <- unique(Xy)
  colnames(uXy) <- pn.p
  uXz <- unique(Xz)
  colnames(uXz) <- pn.phi
  
  ### INITIAL VALUES ###
  if(missing(init.list)){
    beta.z <- rep(0,ncol(Xz))
    beta.y <- rep(0,ncol(Xy))
  }
  else{
    beta.z <- init.list$beta.z
    beta.y <- init.list$beta.y
  }
  
  ###  PRIOR DISTRIBUTIONS ###
  if(missing(prior.list)){
    Q.b.z <- diag(rep(1,ncol(Xz)))
    mu.b.z <- rep(0,ncol(Xz))
    Q.b.y <- diag(rep(1,ncol(Xy)))
    mu.b.y <- rep(0,ncol(Xy))	
  }
  else{
    Q.b.z <- prior.list$phi$Q
    mu.b.z <- prior.list$phi$mu
    Q.b.y <- prior.list$p$Q
    mu.b.y <- prior.list$p$mu
  }
  
  ### STORAGE ###
  beta.z.stor <- matrix(NA, iter, ncol(Xz))
  beta.y.stor <- matrix(NA, iter, ncol(Xy))
  fitted.y.stor <- matrix(NA, iter, nrow(uXy))
  fitted.z.stor <- matrix(NA, iter, nrow(uXz))
  
  colnames(beta.z.stor) <- pn.phi
  colnames(beta.y.stor) <- pn.p
  
  ### BEGIN MCMC ###
  cat("\nprobitCJS MCMC beginning...\n")
  cat("p model = ", as.character(p.model),"\n")
  cat("phi model = ", as.character(phi.model),"\n\n")
  
  tot.iter <- burnin + iter
  st <- Sys.time()
  for(m in 1:tot.iter){
    
    ### UPDATE Z ###
    zvec[idx.z] <- unlist(tapply(1:n.unk.z, id.z, FUN=sample.z, 
                                 mu.y=Xy.z%*%beta.y, mu.z=Xz.z%*%beta.z))
    
    ### UPDATE Z.TILDE ### 
    a <- ifelse(zvec==0, -Inf, 0)
    b <- ifelse(zvec==0, 0, Inf)
    z.tilde <- rtruncnorm(n, a=a, b=b, mean=Xz%*%beta.z, sd=1)
    
    ### BETA.Z UPDATE ###
    idx.z.tilde <- as.logical(ave(zvec, id, FUN=function(x){x + ifelse(diff(c(1,x))==-1, 1, 0)}))
    #idx.z.tilde <- zvec + ifelse(diff(c(1,zvec))==-1 & diff(c(1,id.num))==0, 1, 0)
    V.beta.z.inv <- crossprod(Xz[idx.z.tilde,]) + Q.b.z
    m.beta.z <- solve(V.beta.z.inv, crossprod(Xz[idx.z.tilde,],z.tilde[idx.z.tilde]) + crossprod(Q.b.z,mu.b.z))
    beta.z <- m.beta.z + solve(chol(V.beta.z.inv), rnorm(ncol(Xz),0,1))
    if(m>burnin){
      beta.z.stor[m-burnin,] <- beta.z
      fitted.z.stor[m-burnin,] <- exp(ln.Phi.1(uXz%*%beta.z))
    }
    
    ### UPDATE Y.TILDE ### 
    a <- ifelse(yvec==0, -Inf, 0)
    b <- ifelse(yvec==0, 0, Inf)
    y.tilde <- rtruncnorm(n, a=a, b=b, mean=Xy%*%beta.y, sd=1)
    
    ### BETA.Y UPDATE ###
    V.beta.y.inv <- crossprod(Xy[zvec==1,]) + Q.b.y
    m.beta.y <- solve(V.beta.y.inv, crossprod(Xy[zvec==1,],y.tilde[zvec==1])+crossprod(Q.b.y,mu.b.y))
    beta.y <- m.beta.y + solve(chol(V.beta.y.inv), rnorm(ncol(Xy),0,1))
    if(m>burnin) {
      beta.y.stor[m-burnin,] <- beta.y
      fitted.y.stor[m-burnin,] <- exp(ln.Phi.1(uXy%*%beta.y))
    }
    
    ### TIMING OF SAMPLER ###
    if(m==30){
      tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/30
      ttc <- round((tot.iter-30)*tpi/3600, 2)
      if(ttc>=1) cat("\nApproximate time till completion: ", ttc, " hours\n")
      else cat("\nApproximate time till completion: ", ttc*60, " minutes\n")
    }
    if(100*(m/tot.iter) >= 10 & (100*(m/tot.iter))%%10==0) cat("\n", 100*(m/tot.iter), "% completed\n")
    
  }
  
  ### MAKE SOME OUTPUT ###
  fitted.phi <- mcmc(fitted.z.stor)
  summ.phi <- summary(fitted.phi)
  hpd.phi <- HPDinterval(fitted.phi)
  reals.phi <- data.frame(uXz, 
                          posterior.mean=apply(fitted.z.stor, 2, mean), 
                          CI.lower=hpd.phi[,1], CI.upper=hpd.phi[,2])
  fitted.p <- mcmc(fitted.y.stor)
  summ.p <- summary(fitted.p)
  hpd.p <- HPDinterval(fitted.p)
  reals.p <- data.frame(uXy, 
                        posterior.mean=apply(fitted.y.stor,2,mean), 
                        CI.lower=hpd.p[,1], CI.upper=hpd.p[,2])
  
  
  return(
    list(
      beta.phi=mcmc(beta.z.stor), 
      beta.p=mcmc(beta.y.stor), 
      fitted.p=fitted.p, 
      fitted.phi=fitted.phi,
      reals.phi=reals.phi,
      reals.p=reals.p
      )
    )
}	### END OF FUNCTION ###
