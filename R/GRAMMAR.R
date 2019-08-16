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
#' @param max_itr is specifies the number of permutations.
#' @param num.parallel Number of parallel processes or a predefined socket cluster
#' @param outPath is the file output path 
#' 
#' @return p-value and f-value
#' 
#' @importFrom vegan adonis
#' @importFrom utils read.table
#' @importFrom utils write.table
#' 
#' @examples 
#'    X = as.matrix(read.table("./X_rightdim.txt"))
#'    Y = as.matrix(read.table("./Y_rightdim.txt"))
#'    
#'    K = Kinship(t(X))
#'    VC = varComp(K, Y, X)
#'    
#'    result = run_grammar(K, Y, X, VC, max_itr = 4, num.parallel = 2, outPath = "./")
#'    
#' @export
run_grammar<- function(K, Y, X, VC, max_itr, num.parallel, outPath) {
  ptm <- proc.time()
  
  run_gamma <- function(Y, X, max_itr, num.parallel, outPath) {
    Ng <- dim(X)[2]
    pval <- 1:Ng
    fval <- 1:Ng
    newY <- Y - min(Y)
    
    getp <- function(Y, x, p) {
      res <- vegan::adonis(Y ~ x, perm = p)
      return(res$aov.tab$"Pr(>F)"[1])
    }
    getF <- function(Y, x, p) {
      res <- vegan::adonis(Y ~ x, perm = p)
      return(res$aov.tab$F.Model[1])
    }
    
    esgamma <- function(Y, x, max_itr) {
      Y <- Y
      x <- x
      for (i in 2:max_itr) {
        p <- 10^i
        limit <- 5/p
        pval <- getp(Y, x, p)
        if (pval > limit) {
          break
        }
      }
      return(pval)
    }
    
    cl <- parallel::makeCluster(spec = num.parallel, type = "SOCK", outfile = "log.txt")
    doParallel::registerDoParallel(cl)
    
    '%dopar%' <- foreach::"%dopar%"
    
    
    foreach::foreach(i=1:Ng, .combine = 'rbind', .verbose = T) %dopar% {
      pval[i] <- esgamma(newY, X[, i], max_itr)
      pv <- pval[i]
      fval[i] <- getF(newY, X[, i], 1)
      fv <- fval[i]
      
      cat(i, "\t", pv, "\t", fv, file=paste(outPath, "/result.txt", sep = ""), append = T)
      cat("\n", file=paste(outPath, "/result.txt", sep = ""), sep = "", append=T)
      
    }
    
    tempread <- as.matrix(read.table(paste(outPath, "/result.txt", sep = "")))
    file.remove(paste(outPath, "/result.txt", sep = ""))
    towrite <- tempread[order(tempread[,1]),]
    write.table(towrite, paste(outPath, "/result.txt", sep = ""), row.names = F, col.names = c("Num\t", "P_value\t", "F_value"), quote = F)
    
    parallel::stopCluster(cl)
  }
  
  chol_solve <- function(K) {
    a <- eigen(K)$vectors
    b <- eigen(K)$values
    b[b < 1e-13] <- 1e-13
    b <- 1 / sqrt(b)
    return(a %*% diag(b) %*% t(a))
  }
  
  rotate <- function(Y, sigma) {
    U <- chol_solve(sigma)
    tU <- t(U)
    UY <- tU %*% Y
    return(UY)
  }
  
  snpNum <- dim(X)[2]
  indiNum <- dim(X)[1]
  geneNum <- dim(Y)[2]
  
  Vg <- median(VC[,1])		# Variance components
  Ve <- median(VC[,2])
  
  I <- diag(indiNum)
  sigma <- Vg*K + Ve*I
  
  UY <- rotate(Y,sigma)		# Rotate genotypes and phenotypes
  UX <- rotate(X,sigma)
  
  pf <- run_gamma(UY, UX, max_itr, num.parallel, outPath)
  
  calTime <- proc.time() - ptm
  write.table(calTime, paste(outPath, "/time.txt", sep = ""), row.names = F, col.names = F, quote = F)
}