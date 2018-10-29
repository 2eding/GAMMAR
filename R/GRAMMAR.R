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
#' run_gamma function calculates p-value and f-value.
#' 
#' @param K is obtained from the Kinship function
#' @param UY is obtained from the varComp function
#' @param UX is obtained from the varComp function
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
#'    ps = run_gamma(Y, X, max_itr = 4, num.parallel = 2)
#'    # The default value for num.parallel is 2
#'    
#'    # ps[1] = p-value
#'    # ps[2] = f-value
#' @export
run_gamma<- function(UY, UX, max_itr = 4, num.parallel = 2) {
  
  getp <- function(Y, x, p, num.parallel = 2) {
    res = vegan::adonis(Y ~ x, perm = p, parallel = num.parallel)
    return(res$aov.tab$"Pr(>F)"[1])
  }
  getF <- function(Y, x, p, num.parallel = 2) {
    res = vegan::adonis(Y ~ x, perm = p, parallel = num.parallel)
    return(res$aov.tab$F.Model[1])
  }

  gamma <- function(Y, x, max_itr = 4, num.parallel = 2) {
    for (i in 2:max_itr) {
      p = 10^i
      limit = 5/p
      pval = getp(Y, x, p, num.parallel)
      if (pval > limit) {
        break
      }
    }
    return(pval)
  }
  
  Ng = dim(X)[2]
  pval = 1:Ng
  fval = 1:Ng
  newY = Y-min(Y)
  for (i in 1:Ng) {
    pval[i] = gamma(newY, X[,i], max_itr, num.parallel)
    fval[i] = getF(newY,X[,i],1, num.parallel)
    cat(i,". f=",fval[i]," p=",pval[i],"\n")
  } 
  write.table(pval, "P.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(fval, "F.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(list(pval, fval))
}



