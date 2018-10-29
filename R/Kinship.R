#' Kinship
#'
#' Kinship function will estimate kinship coefficient
#' 
#' @param X is the SNP matrix, individual x snp
#' 
#' @examples
#'    X = as.matrix(fread("SNP_rightdim.txt"))
#'    K = Kinship(t(X))
#' @export

Kinship <- function(X){
  X_norm = ((X - mean(X, na.rm = T)) / sd(X, na.rm = T))
  return (cor(X_norm))
}