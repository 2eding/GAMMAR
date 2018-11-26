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
#' 
#' @examples 
#'    X = as.matrix(fread("X_rightdim.txt"))
#'    K = Kinship(t(X))
#'    Y = as.matrix(fread("Y_rightdim.txt"))
#'    VC = varComp(K, Y, X)
#' @export


varComp <- function(K, Y, X){
  ptm <- proc.time()
  snpNum <- dim(X)[2]
  indiNum <- dim(X)[1]
  geneNum <- dim(Y)[2]
  K_origin = K
  
  for (i in 1:dim(Y)[2]) {
    K = K_origin
    print(i)
    
    e = lmmlite::eigen_rotation(K, t(Y)[i,], use_cpp = T)
    VC = lmmlite::fitLMM(
      e$Kva,
      e$y,
      e$X,
      reml = T,
      use_cpp = T,
      tol = 1e-6,
      check_boundary = T
    )
    
    write.table(
      VC$sigmasq_g,
      "Vg_temp.txt",
      row.names = F,
      col.names = F,
      append = T,
      quote = F,
      sep = "\n"
    )
    write.table(
      VC$sigmasq_e,
      "Ve_temp.txt",
      row.names = F,
      col.names = F,
      append = T,
      quote = F,
      sep = "\n"
    )
  }
  
  VCbind = rbind("Vg.txt", "Ve.txt")
  file.remove("Vg.txt")
  file.remove("Ve.txt")
  write.table(VCbind, "VC.txt", row.names = F, col.names = F, quote = F)
  Vg = median(VCbind[1])
  Ve = median(VCbind[2])
  
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
  
  print(proc.time() - ptm)
  return(list(
    "Vg" = Vg,
    "Ve" = Ve,
    "UY" = UY,
    "UX" = UX
  ))
}
