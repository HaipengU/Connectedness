#' Computation of genomic relationship matrix with SNP markers.
#'
#' Use genomic relationship matrix to compute genomic connectedness.
#' 
#' @param Wmatrix a raw SNP marker matrix with dimension n by m, where n is individual and m is marker. 
#' @param maf a minor allele frequency used for quality control; e.g., 0.05
#' @param type a type of genomic relationship matrix 
#' @return a n by n genomic relationship matrix.
#' 
#' @examples 
#' computeG()
#' 
#' @export
#' 
computeG <- function(Wmatrix, maf, type) {
  p1 <- (colMeans(Wmatrix, na.rm = T)/2)
  set.seed(0)
  for (j in 1:ncol(Wmatrix)) {
    Wmatrix[,j] <- ifelse(is.na(Wmatrix[,j]), rbinom(1,2,p1[j]), Wmatrix[,j])
  }
  p2 <- colMeans(Wmatrix) / 2
  maf2 <- pmin(p2, 1-p2)
  maf.index <- which(maf2 < maf)
  W <- Wmatrix[, -maf.index]
  if (type == 'G1') { 
    W_c <- scale(W, center = TRUE, scale = FALSE)
    G1 <- tcrossprod(W_c) / sum(2 * p2 * (1 - p2))
    return(G1)
  } else if (type == 'G2'){
    W_cs <- scale(W)
    G2 <- tcrossprod(W_cs) / ncol(W_cs)
    return(G2)
  } else if (type == 'G0.5') {
    W_cs <- (W - 2 * 0.5) / sqrt(2 * 0.5 * (1 - 0.5)) 
    G0.5 <- tcrossprod(W_cs) / ncol(W_cs)
    return(G0.5)
  } else if (type == 'G1_s') {
    W_c <- scale(W, center = TRUE, scale = FALSE)
    G1 <- tcrossprod(W_c) / sum(2 * p2 * (1 - p2))
    G1_s<- 2 * (G1 - min(G1)) / (max(G1) - min(G1))
    return(G1_s)
  } else if (type == 'G2_s') {
    W_cs <- scale(W)
    G2 <- tcrossprod(W_cs) / ncol(W_cs)
    G2_s<- 2 * (G2 - min(G2)) / (max(G2) - min(G2))
    return(G2_s)
  }  
} 

