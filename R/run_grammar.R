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
#' @param outPath is a parameter that specifies the path of the result file
#' @param name is a parameter that specifies the name of the result file
#' 
#' @return p-value and f-value
#' 
#' @importFrom vegan adonis
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' 
#' @examples 
#'    X <- as.matrix(data.table::fread(Genotypesdata))
#'    Y <- as.matrix(data.table::fread(Phenotypesdata))
#'    
#'    K <- Kinship(X)
#'    VC <- varComp(K, Y, X)
#'    
#'    result <- run_grammar(K, Y, X, VC, max_itr = 4, num.parallel = 2, outPath = "./testdir/", name = "result.txt")
#'    
#' @export
run_grammar<- function(K, Y, X, VC, max_itr, num.parallel, outPath, name) {
  ptm <- proc.time()
  
  run_gamma <- function(Y, X, max_itr, num.parallel, outPath, name) {
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
    
    cl <- parallel::makeCluster(spec = num.parallel, type = "SOCK")
    doParallel::registerDoParallel(cl)
    
    '%dopar%' <- foreach::"%dopar%"
    
    
    foreach::foreach(i=1:Ng, .combine = 'rbind', .verbose = T) %dopar% {
      pval[i] <- esgamma(newY, X[, i], max_itr)
      pv <- pval[i]
      fval[i] <- getF(newY, X[, i], 1)
      fv <- fval[i]
      
      saveresult <- c(i, "\t", pv, "\t", fv, "\n")
      Sys.sleep(0.0001)
      cat(saveresult, file=paste(outPath, "/", name, sep = ""), append=T)
    }
    
    tempread <- as.matrix(read.table(paste(outPath, "/", name, sep = "")))
    file.remove(paste(outPath, "/", name, sep = ""))
    towrite <- tempread[order(tempread[,1]),]
    write.table(towrite, paste(outPath, "/", name, sep = ""), row.names = F, col.names = c("SNP_Num\t", "P_value\t", "F_value"), quote = F)
    
    parallel::stopCluster(cl)
    
    return(towrite)
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
  
  grammar_result <- run_gamma(UY, UX, max_itr, num.parallel, outPath, name)

  print(proc.time() - ptm)
  return(grammar_result)
}