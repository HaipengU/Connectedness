var_est <- function(y,K,X){
  library(rrBLUP)
  GBLUP <- mixed.solve(y = y, K = K, X = X)
  varcomp <- list(sigma2a = GBLUP$Vu, sigma2e = GBLUP$Ve)
  return(varcomp)
}
