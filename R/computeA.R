#' Computation of numerator relationship matrix.
#'
#' Use sorted pedigree to calculate numerator relationship matrix. The missing parents are coded as 0. 
#' 
#' @param Progeny a vector of sorted animals
#' @param Sire a vector of Sires according to progenies, where the missing Sire is coded as 0
#' @param Dam a vector of Dams according to progenies, where the missing Dam is coded as 0
#' @return a n by n numerator relationship matrix.
#' 
#' @examples 
#' computeA()
#' 
#' @export
#' 
computeA <- function(Progeny, Sire, Dam){
  if (any(duplicated(Progeny))) stop("Progeny must be unique")
  n <- length(Progeny)
  A <- diag(n)
  for (i in 1 : n){
    if (Sire[i] == 0 && Dam[i] != 0){
      for (j in 1 : i-1){
        A[i, j] <- A[j, i] <- 0.5 *(A[j, Dam[i]]) 
      } 
    } else if (Sire[i] != 0 && Dam[i] == 0) {
      for (j in 1 : i-1){
        A[i, j] <- A[j, i] <- 0.5 *(A[j, Sire[i]])
      }
    } else if (Sire[i] != 0 && Dam[i] != 0){
      for (j in 1 : i -1){
        A[i, j] <- A[j, i] <- 0.5 *(A[j, Sire[i]] + A[j, Dam[i]])
      }
      A[i, i] <- A[i, i] + 0.5 * (A[Sire[i], Dam[i]])
    }
  }
  return(A)
}


