#' Estimates of variance componets. 
#'
#' Maximum Likelihood estimation of variance components with eigenvectors.
#' 
#' @param y a n by 1 vector of phenotypes. 
#' @param Evector a n by n matrix with columns according to eigenvectors.
#' @param Evalue a vector contains the n eigenvalues.
#' 
#' @return A list of variance components.
#' 
#' @examples 
#' varcomp()
#' 
#' @export
#' 
varcomp <- function(y, Evector, Evalue){
  startVal <- c(0.5, 0.5)
  var.opt <- optim(fn = log.Lik, y=y, Evector = Evector, Evalue = Evalue, par = startVal,
                   hessian=FALSE) 
  var.est <- list(Ve = var.opt$par[1] , Vu = var.opt$par[2])
  return(var.est)
}


# loglikelihood func
log.Lik <- function(y, Evector, Evalue, startVar){
  y <- y - mean(y)
  n <- length(y)
  V <- Evector
  D <- Evalue
  V_y <- crossprod(V,y) 
  V_2y <- as.vector(V_y)^2 
  sigma2e <- startVar[1] 
  sigma2a <- startVar[2] 
  lambda <- sigma2a/sigma2e 
  DStar <- (D * lambda + 1)
  sumLogD <- log(prod(DStar))
  part1 <- ( n * log(sigma2e) + sumLogD ) 
  part2 <- (sum(V_2y/DStar)) / sigma2e
  LogLik <- part1 + part2
  return(LogLik)
}
