#' Create release-recapture history data for MCMC
#' 
#' Creates needed constructs from the release-recapture history.
#' 
#' 
#' @param data A data frame where each row contains the capture history and covariates for each individual. The cature history must be spread over multiple columns
#' @param id.col Integer or character indicating the column that holds a personal ID for each individual
#' @param ch.cols Vector of integers giving the columns that contain the capture histories
#' @param covariate.cols Vector of integers that indicate the columns containg covariate information
#' @return A data frame 
#' @export
#' @author Devin Johnson

make.probitCJS.data <- function(data, id.col, ch.cols, covariate.cols){
  if(missing(id.col)){
    data$id <- factor(1:nrow(data))
    id.col <- ncol(data)
  }
  data <- data[order(data[,id.col]),]
  id <- data[,id.col]
  ch <- as.matrix(data[,ch.cols])
  if(any(!ch%in%c(0,1))) stop("There are values other than 0 and 1 in the capture history columns\n")
  data.out <- data[,c(id.col,covariate.cols)]
  f = ncol(ch)-apply(col(ch)[,ncol(ch):1]*ch,1,max)+1
  l = apply(col(ch)*ch,1,max)
  ch[col(ch)<=f] <- NA
  num.smp <- ncol(ch)-f
  Z <- as.matrix(sweep(col(ch),1,f,'>=')*sweep(col(ch),1,l,'<='))
  Z <- ifelse(sweep(col(ch),1,l,'>'),NA,Z)
  YZ <- cbind(Y=as.vector(t(ch)),Z=as.vector(t(Z)))
  idx <- !is.na(YZ[,1])
  YZ <- YZ[idx,]
  Cohort <- f
  cohort <- factor(f)
  data.out <- cbind(data.out, cohort, Cohort)[rep(1:length(id),num.smp),]
  Time <- as.vector(t(col(ch)))[idx]-2
  time <- factor(Time+2)
  Age <- as.vector(t(col(ch) - f))[idx]
  age <- factor(Age)
  data.out <- cbind(data.out, age, Age, time, Time, YZ)
  attr(data.out, "id") <- names(data)[id.col]
  class(data.out) <- c(class(data.out),"probitCJS.data")
  return(data.out)
}
