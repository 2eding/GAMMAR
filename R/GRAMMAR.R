#######################################################################################################################
#######################################################################################################################
#################################### To use RpackGAMMA, need to R version >= 3.5.0 ####################################
#################### If you have any question using RpackGAMMA, please email me 2eding@dongguk.edu ####################
#######################################################################################################################
#######################################################################################################################
#################### input X, individual x snp ################### input Y, individual x phenotype ####################
###################################### Check installed.package to use RpackGAMMA ######################################


################################################## Run GAMMA ###################################################

#' Generalized Renown Analysis of Molecular variance for Mixed model Analysis R package
#' 
#' run_grammar function calculates p-value and f-value.
#' 
#' @param K is obtained from the Kinship function
#' @param Y is the phenotype matrix, individual x phenotype
#' @param X is the SNP matrix, individual x snp
#' @param VC is obtained from the varComp function
#' @param mat_itr is specifies the number of permutations.
#' @param num.parallel Number of parallel processes or a predefined socket cluster
#' 
#' @return p-value and f-value
#' 
#' @importFrom 'vegan::' code before adonis
#' @examples 
#'    X = as.matrix(fread("SNP_rightdim.txt"))
#'    K = Kinship(t(X))
#'    
#'    Y = as.matrix(fread("Phenotype_rightdim.txt"))
#'    VC = varComp(K, Y, X)
#'    
#'    ps = run_grammar(VC$UY, VC$UX, max_itr = 4, num.parallel = 2)
#'    
#'    # ps[1] = p-value
#'    # ps[2] = f-value
#' @export
run_grammar<- function(K, Y, X, VC, max_itr, num.parallel) {
  ptm <- proc.time()

  getp <- function(Y, x, p) {
    require(parallel)
    res = vegan::adonis(Y ~ x, perm = p, parallel = getOption("mc.cores"))
    return(res$aov.tab$"Pr(>F)"[1])
  }
  getF <- function(Y, x, p) {
    require(parallel)
    res = vegan::adonis(Y ~ x, perm = p, parallel = getOption("mc.cores"))
    return(res$aov.tab$F.Model[1])
  }

  gamma <- function(Y, x, max_itr) {
    Y = Y
    x = x
    for (i in 2:max_itr) {
      p = 10^i
      limit = 5/p
      pval = getp(Y, x, p)
      if (pval > limit) {
        break
      }
    }
    return(pval)
  }

  run_gamma <- function(Y, X, max_itr, num.parallel) {
    Ng = dim(X)[2]
    pval = 1:Ng
    fval = 1:Ng
    newY = Y - min(Y)
    
    require(snow)
    cl <- parallel::makeCluster(spec = num.parallel, type = "SOCK")
    doParallel::registerDoParallel(cl)
    Sys.setenv("MC_CORES"=num.parallel)
    
    '%dopar%' <- foreach::"%dopar%"
    
    foreach::foreach(i=1:Ng) %dopar% {
      pval[i] = gamma(newY, X[, i], max_itr)
      fval[i] = getF(newY, X[, i], 1)
      cat(i, ". f =", fval[i], " p =", pval[i], "\n")
      write.table(pval[i], "P.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, append=T, sep="\n")
      write.table(fval[i], "F.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, append=T, sep="\n")
    }
    
    # for (i in 1:Ng) {
    #   require(parallel)
    #   pval[i] = gamma(newY, X[, i], max_itr)
    #   fval[i] = getF(newY, X[, i], 1)
    #   cat(i, ". f =", fval[i], " p =", pval[i], "\n")
    #   write.table(pval[i], "P.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, append=T, sep="\n")
    #   write.table(fval[i], "F.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, append=T, sep="\n")
    # }
    # return(list("p" = pval, "f" = fval))
  }
  
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
  
  snpNum <- dim(X)[2]
  indiNum <- dim(X)[1]
  geneNum <- dim(Y)[2]
  
  Vg = median(VC[,1])		# Variance components
  Ve = median(VC[,2])
  
  I = diag(indiNum)
  sigma = Vg*K + Ve*I
  
  UY = rotate(Y,sigma)		# Rotate genotypes and phenotypes
  UX = rotate(X,sigma)
  
  pf = run_gamma(UY, UX, max_itr, num.parallel)

  parallel::stopCluster(cl)
  
  print(proc.time() - ptm)
  return(list("P" = pf$p, "F" = pf$f))
}