#' Estimates of variance componets. 
#'
#' Maximum Likelihood estimation of variance components using the eigenvalues and eigenvectors, 
#'   which are derived from the eigendecomposition of the relationship matrix. 
#' 
#' @param y A vector includes the phenotypes of n individuals.  
#' @param Evector A matrix  (n x n) with columns according to eigenvectors from the eigendecomposition of the relationship matrix. 
#' @param Evalue A vector which contains the eigenvalues of n individuals from the eigendecomposition of the relationship matrix. 
#' 
#' @return A list contains variance components
#' \describe{
#'   \item{$Ve}{An estimate of residual variance.}
#'   \item{$Vu}{An estimate of additive genetic variance.}
#' }
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' 
#' @example man/examples/varcomp.R
#' 
#' @export
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
