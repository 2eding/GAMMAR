#' varComp
#'
#' varComp function will estimate variance components
#' 
#' @param K is obtained from the Kinship function
#' @param Y is the phenotype matrix, individual x phenotypes
#' @param num.parallel Number of parallel processes or a predefined socket cluster
#' 
#' @return Variance components Vg: genetic factor Ve: environment factor
#' 
#' @importFrom pbmcapply pbmclapply
#' @importFrom regress regress
#' 
#' @examples 
#'    X <- as.matrix(data.table::fread(Genotypesdata))
#'    Y <- as.matrix(data.table::fread(Phenotypesdata))
#'    K <- Kinship(X)
#'    
#'    VC <- varComp(Y, K, num.parallel)
#' @export


varComp <- function(Y, K, num.parallel){
  ptm <- proc.time()

  Y <- t(Y)
  calvc <- function(Y, K, core){
    list_Y <- as.list(as.data.frame(t(Y)))
    length(list_Y)
    
    raw_vc <- pbmcapply::pbmclapply(list_Y, function(Y){
      regress::regress(Y ~ 1, ~ K, pos=c(T, T), tol=1e-8)$sigma},
      mc.cores = core
      )
    
    return(do.call(rbind, raw_vc))
  }
  
  ncore <- num.parallel
  vc <- calvc(Y, K, core = ncore)
  write.table(vc, "./VC.txt", row.names = F, col.names = F, quote = F)
  
  print(proc.time() - ptm)
  return(vc)
}