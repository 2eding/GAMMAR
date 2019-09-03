#' Kinship
#'
#' Kinship function will estimate kinship coefficient
#' 
#' @param X is the SNP matrix, individual x SNPs
#' 
#' @importFrom stats sd
#' @importFrom stats cor
#' @importFrom utils write.table
#' 
#' @examples
#'    X = as.matrix(read.table(Genotypesdata))
#'    K = Kinship(t(X))
#' @export

Kinship <- function(X){
  ptm <- proc.time()

  varjj <- function(data){
    return(sum((mean(data) - data)^2)/length(data))
  }
  
  W <- t(X)
  n <- nrow(W)
  m <- ncol(W)
  keep <- c()
  
  # i=1
  for (i in 1:m){
    mn <- mean(W[!is.na(W[,i]),i])
    W[is.na(W[,i]),i] <- mn
    vr <- varjj(W[,i])
    
    if (vr == 0) s <- 1.0
    else s <- sqrt(vr)
    
    keep <- c(keep,i)
    W[,i] <- (W[,i] - mn) / s
  }
  W <- W[, keep]
  # K <- (W %*% t(W)) * 1.0/m
  K <- tcrossprod(W) * 1.0/m
  write.table(K, "./K.txt", row.names = F, col.names = F, quote = F)
  print(proc.time() - ptm)
  return(K)
}