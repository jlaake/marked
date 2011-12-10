#' Split a character string of capture histories into a matrix
#'
#' Split a character string of capture histories into a matrix. The matrix is appended to the original data set after
#'
#' @param x A vector containing the character strings of capture histories or the column number or name in the data set \code{data}
#' @param data A data fram containing x if x indicates a column in a data frame
#' @return A data frame
#' @export
#' @author Devin Johnson

splitCH <- function(x="ch", data=NULL){
  if((is.character(x) & length(x)==1) & !is.null(data)){
    ch <- as.character(data[,x])
    chmat <- t(sapply(ch, FUN=function(v){as.numeric(unlist(strsplit(v, split="")))}))
    colnames(chmat) <- paste("Time", c(1:ncol(chmat)),sep="")
    rownames(chmat) <- NULL
    return(cbind(data,chmat))
  }
  else{
    chmat <- t(sapply(x, FUN=function(v){as.numeric(unlist(strsplit(v, split="")))}))
    colnames(chmat) <- paste("Time", c(1:ncol(chmat)),sep="")
    return(chmat)
  }
}
