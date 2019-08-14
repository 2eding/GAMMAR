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
#'    X = as.matrix(read.table("X_rightdim.txt"))
#'    K = Kinship(t(X))
#' @export

Kinship <- function(X){
  ptm <- proc.time()
  X_norm = (X - mean(X, na.rm = T) / sd(X, na.rm = T)) # normalize
  kin = cor(X_norm)
  write.table(kin, "K.txt", row.names = F, col.names = F, quote = F)
  print(proc.time() - ptm)
  return (kin)
}