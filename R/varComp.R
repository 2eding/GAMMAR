#' varComp
#'
#' varComp function will estimate variance components
#' 
#' @param K is obtained from the Kinship function
#' @param Y is the phenotype matrix, individual x phenotype
#' @param X is the SNP matrix, individual x snp
#' 
#' @return Variance components Vg: genetic factor(VC[1]) Ve: environment factor(VC[2])
#' @return UY: rotate(Y, sigma)(VC[3]) UX: rotate(X, sigma)(VC[4])
#' 
#' @importFrom 'lmmlite::' code before eigen_rotation
#' @importFrom 'lmmlite::' code before fitLMM
#' @importFrom 'progress::' code before progress_bar
#' @importFrom 'data.table::' code before fread
#' 
#' @examples 
#'    X = as.matrix(fread("SNP_rightdim.txt"))
#'    K = Kinship(t(X))
#'    Y = as.matrix(fread("Phenotype_rightdim.txt"))
#'    VC = varComp(K, Y, X)
#' @export


varComp <- function(K, Y, X){
  pvc = progress::progress_bar$new(format=" [:bar] :percent", total = dim(Y)[2]*dim(X)[2])
  for(i in 1:dim(Y)[2]) {
    for(j in 1:dim(X)[2]){
      e = lmmlite::eigen_rotation(K, t(Y)[i,], X[,j], use_cpp = T)
      VC = lmmlite::fitLMM(e$Kva, e$y, e$X, reml = T,use_cpp = T, tol = 1e-6, check_boundary = T)
      pvc$tick()
      write.table(VC$sigmasq_g, "sigmasq_g1.txt", row.names = F, col.names = F, append = T, quote = F, sep="\n")
      write.table(VC$sigmasq_e, "sigmasq_e1.txt", row.names = F, col.names = F, append = T, quote = F, sep="\n")
    }
    Vg_temp = median(as.matrix(data.table::fread("sigmasq_g1.txt")))
    Ve_temp = median(as.matrix(data.table::fread("sigmasq_e1.txt")))
    write.table(Vg_temp, "Vg_temp.txt", row.names = F, col.names = F, quote = F, append = T, sep="\n")
    write.table(Ve_temp, "Ve_temp.txt", row.names = F, col.names = F, quote = F, append = T, sep="\n")
    file.remove("sigmasq_g1.txt")
    file.remove("sigmasq_e1.txt")
  }
  Vg = median(as.matrix(data.table::fread("Vg_temp.txt")))
  Ve = median(as.matrix(data.table::fread("Ve_temp.txt")))

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
  
  return(list(Vg, Ve, UY, UX))
}
