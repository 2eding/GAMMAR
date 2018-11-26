#' varComp
#'
#' varComp function will estimate variance components
#' 
#' @param K is obtained from the Kinship function
#' @param Y is the phenotype matrix, individual x phenotype
#' @param X is the SNP matrix, individual x snp
#' 
#' @return Variance components Vg: genetic factor(VC[1]) Ve: environment factor(VC[2])
#' @return UY is the third output value of the varComp function.
#' @return UX is the fourth output value of the varComp function.
#' 
#' @importFrom 'lmmlite::' code before eigen_rotation
#' @importFrom 'lmmlite::' code before fitLMM
#' @importFrom 'progress::' code before progress_bar
#' @importFrom 'data.table::' code before fread
#' 
#' @examples 
#'    X = as.matrix(fread("X_rightdim.txt"))
#'    K = Kinship(t(X))
#'    Y = as.matrix(fread("Y_rightdim.txt"))
#'    VC = varComp(K, Y, X)
#' @export


varComp <- function(K, Y, X){
  snpNum <- dim(X)[2]
  indiNum <- dim(X)[1]
  geneNum <- dim(Y)[2]
  X0 = array(1, c(indiNum, geneNum))
  X0_origin = X0
  K_origin = K
  
  for (i in 1:dim(Y)[2]) {
    ptm <- proc.time()
    pvc = progress::progress_bar$new(format = " [:bar] :percent", total = geneNum)
    
    X0 = X0_origin
    K = K_origin
    
    e = lmmlite::eigen_rotation(K, t(Y)[i,], X0, use_cpp = T)
    VC = lmmlite::fitLMM(
      e$Kva,
      e$y,
      e$X,
      reml = T,
      use_cpp = T,
      tol = 1e-6,
      check_boundary = T
    )
    pvc$tick()
    
    write.table(
      VC$sigmasq_g,
      "Vg.txt",
      row.names = F,
      col.names = F,
      append = T,
      quote = F,
      sep = "\n"
    )
    write.table(
      VC$sigmasq_e,
      "Ve.txt",
      row.names = F,
      col.names = F,
      append = T,
      quote = F,
      sep = "\n"
    )
  }
  
  Vg = median(as.matrix(fread("Vg.txt")))
  Ve = median(as.matrix(fread("Ve.txt")))
  
  chol_solve <- function(K) {
    a = eigen(K)$vectors
    b = eigen(K)$values
    b[b < 1e-13] = 1e-13
    b = 1 / sqrt(b)
    return(a %*% diag(b) %*% t(a))
  }
  rotate <- function(Y, sigma) {
    U <- chol_solve(sigma)
    tU <- t(U)
    UY = tU %*% Y
    return(UY)
  }
  
  I = diag(indiNum)
  sigma = Vg * K + Ve * I
  UY = rotate(Y, sigma)		# Rotate genotypes and phenotypes
  UX = rotate(X, sigma)
  
  Sys.sleep(1 / 100)
  pvc$terminate()
  
  print(proc.time() - ptm)
  return(list(
    "Vg" = Vg,
    "Ve" = Ve,
    "UY" = UY,
    "UX" = UX
  ))
}
