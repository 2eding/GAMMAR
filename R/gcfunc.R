#' gcfunc
#' 
#' Perform genomic control of p value
#' 
#' @param P is the p value from the result of GRAMMAR function
#' 
#' @return completed value of the genomic control of p value
#' 
#' @importFrom stats median
#' @importFrom stats pchisq
#' @importFrom stats qchisq
#' 
#' @examples 
#'    P = as.matrix(fread("P.txt"))
#'    gcp = gcfunc(P)
#' @export
gcfunc <- function(P){
  ps = as.matrix(P)
  lambda = median(qchisq(ps,1, lower.tail=F),na.rm=TRUE) / qchisq(0.5,1)
  gc <- function(ps, lambda) {
    return(pchisq(qchisq(ps, 1, lower.tail=F)/lambda, 1, lower.tail=F))
  }
  p = gc(ps, lambda)
  
  return(p)
}