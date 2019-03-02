#' Estimate variance componets. 
#'
#' Use rrBLUP to estimate variance components. 
#' 
#' @param y a n by 1 vector of phenotypes. 
#' @param K a relationship matrix with a dimension of n by n.
#' @param X a design matrix which associates fixed effects with y.
#' 
#' @return a lit of variance components.
#' 
#' @examples 
#' VarComp()
#' 
#' @export
#' 
VarComp <- function(y,K,X){
  library(rrBLUP)
  GBLUP <- mixed.solve(y = y, K = K, X = X)
  varcomp <- list(sigma2a = GBLUP$Vu, sigma2e = GBLUP$Ve)
  return(varcomp)
}
