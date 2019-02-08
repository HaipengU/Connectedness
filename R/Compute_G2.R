#' Computation of genomic relationship matrix with SNP markers.
#'
#' Use genomic relationship matrix to compute genomic connectedness.
#' 
#' @param Wmatrix a raw SNP marker matrix with dimension n by m, where n is individual and m is marker. 
#' @param maf a minor allele frequency used for quality control; e.g., 0.05
#' 
#' @return a n by n genomic relationship matrix.
#' 
#' @examples 
#' computeG2()
#' 
#' @export
#' 
computeG2 <- function(Wmatrix, maf) {
    set.seed(1213)
    p1 <- (colMeans(Wmatrix, na.rm = T)/2)
    # imputation for missing markers with rbinom
    for (j in 1:ncol(Wmatrix)){
      Wmatrix[,j] <- ifelse(is.na(Wmatrix[,j]), rbinom(1,2,p1[j]), Wmatrix[,j])
    }
    # remove snp with less than 0.05 MAF
    p2 <- colMeans(Wmatrix)/2
    maf2 <- pmin(p2,1-p2)
    maf.index <- which(maf2 < maf)
    W <- Wmatrix[,-maf.index]#
    W_cs <- scale(W)
    Gmatrix <- tcrossprod(W_cs)/ncol(W_cs)
    return(Gmatrix)
} 
