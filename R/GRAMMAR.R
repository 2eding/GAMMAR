#' Generalized Renown Analysis of Molecular variance for Mixed model Analysis R package
#' 
#' GRAMMAR function calculates p-value and f-value just one step
#' 
#' @param X is the SNP matrix, individual x snp
#' @param Y is the phenotype matrix, individual x phenotype
#' @param max_itr is specifies the number of permutations.
#' @param num.parallel Number of parallel processes or a predefined socket cluster
#' @param outPath is the result file output path 
#' 
#' @return p-value and f-value
#' 
#' @importFrom data.table fread
#' @importFrom utils write.table
#' 
#' @examples 
#'    
#'    R CMD BATCH --args -Xpath= -Ypath= -max_itr= -num.parallel= -outPath=./ GRAMMAR.R GRAMMAR.log
#'    
#' @export


GRAMMAR <- function(X, Y, max_itr, num.parallel, outPath) {
  ptm <- proc.time()
  
  args <- commandArgs(TRUE)
  
  if(length(args) != 5){
    print("Usage: R CMD BATCH --args -Xpath= -Ypath= -max_itr= -num.parallel= -outPath= GRAMMAR.R GRAMMAR.log")
    stop()
  }
  
  for(i in 1:5){
    paramName = strsplit(args[i], "=")[[1]][1]
    param = strsplit(args[i], "=")[[1]][2]
    if(paramName == "-Xpath")
      Xpath = param
    else if(paramName == "-Ypath")
      Ypath = param
    else if(paramName == "-max_itr")
      max_itr = param
    else if(paramName == "-num.parallel")
      num.parallel = param
    else if(paramName == "-outPath")
      outPath = param
    else{
      cat("Error: Wrong parameter name ", paramName)
      stop()
    }
  }
  
  X <- as.matrix(data.table::fread(Xpath))
  Y <- as.matrix(data.table::fread(Ypath))
  K <- GRAMMAR::Kinship(t(X))
  VC <- GRAMMAR::varComp(K, Y, X)
  
  Result <- GRAMMAR::run_grammar(K, Y, X, VC, max_itr, num.parallel, outPath)
  
  write.table(Result[,6], paste(outPath, "/P.txt", sep=""), row.names=F, col.names=F, quote=F)
  write.table(Result[,4], paste(outPath, "/F.txt", sep=""), row.names=F, col.names=F, quote=F)
  
  cat(proc.time() - ptm)
}