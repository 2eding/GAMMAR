# 
# if("devtools" %in% rownames(installed.packages()) == FALSE) {
#   install.packages("devtools")
# }
# if("data.table" %in% rownames(installed.packages()) == FALSE) {
#   install.packages("data.table")
# }
# if("lmmlite" %in% rownames(installed.packages()) == FALSE){
#   install_github("kbroman/lmmlite")
# }
# if("RcppEigen" %in% rownames(installed.packages()) == FALSE){
#   install.packages("RcppEigen")
# }
# if("Rcpp" %in% rownames(installed.packages()) == FALSE){
#   install.packages("Rcpp")
# }
# if("progress" %in% rownames(installed.packages()) == FALSE){
#   install.packages("progress")
# }
# 
# library(devtools);
# library(data.table);
# library(lmmlite);
# library(RcppEigen);
# library(Rcpp);
# library(progress);

#' varComp
#'
#' varComp function will estimate variance components
#' 
#' @param K is obtained from the Kinship function
#' @param Y is the phenotype matrix, individual x phenotype
#' @param X is the SNP matrix, individual x snp
#' 
#' @return Variance components Vg: genetic factor Ve: environment factor
#' 
#' @importFrom 'lmmlite::' code before lmmlite
#' 
#' @examples 
#'    X = as.matrix(fread("SNP_rightdim.txt"))
#'    K = Kinship(t(X))
#'    Y = as.matrix(fread("Phenotype_rightdim.txt"))
#'    VC = varComp(K, Y, X)
#' @export


varComp <- function(K, Y, X){
  pvc = progress_bar$new(format=" [:bar] :percent", total = dim(Y)[2]*dim(X)[2])
  for(i in 1:dim(Y)[2]) {
    for(j in 1:dim(X)[2]){
      e = lmmlite::eigen_rotation(K, t(Y)[i,], X[,j], use_cpp = T)
      VC = lmmlite::fitLMM(e$Kva, e$y, e$X, reml = T,use_cpp = T, tol = 1e-6, check_boundary = T)
      pvc$tick()
      write.table(VC$sigmasq_g, "sigmasq_g1.txt", row.names = F, col.names = F, append = T, quote = F, sep="\n")
      write.table(VC$sigmasq_e, "sigmasq_e1.txt", row.names = F, col.names = F, append = T, quote = F, sep="\n")
    }
    Vg = median(as.matrix(fread("sigmasq_g1.txt")))
    Ve = median(as.matrix(fread("sigmasq_e1.txt")))
    # file.remove("sigmasq_g1.txt")
    # file.remove("sigmasq_e1.txt")
    list(Vg, Ve)
  }
  
  chol_solve <- function(K) {
    a = eigen(K)$vectors
    b = eigen(K)$values
    b[b<1e-13] = 1e-13
    b = 1/sqrt(b)
    return(a%*%diag(b)%*%t(a))
  }
  rotate <- function(Y, sigma) {
    U <- chol_solve(sigma)
    tU <-t(U)
    UY = tU%*%Y
    return(UY)
  }
  
  Sys.sleep(1/100)
  pvc$terminate()
  snpNum <- dim(X)[2]
  indiNum <- dim(X)[1]
  geneNum <- dim(Y)[2]
  I = diag(indiNum)
  sigma = Vg*K + Ve*I
  UY = rotate(Y, sigma)
  UX = rotate(X, sigma)
}
