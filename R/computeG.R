#' Compute genomic relationship matrix
#'
#' Use SNP markers to derive additive genomic relationship matrix. Missing marker is allowed and should be coded as NA. 
#' 
#' @param snpmatrix A marker matrix with a dimension of n by m and the elements are coded as 0, 1, 2 or NA, 
#'   where n and m indicate the total number of individuals and markers, accordingly.  
#' @param maf A minor allele frequency for quality control (e.g., 0.05).
#' @param impute Imputation method for missing markers if applicable. Two methods of 'mean' and 'rbinom' are available, 
#'   where the  'mean' imputes the missing marker using mean, and 'rbinom' imputes the missing marker by radom sampling from 
#'   a binomial distribution. The default method is 'mean'. This argument will be ignored the \code{snpmatrix} does not inlcude missing markers. 
#' @param method A type of genomic relationship matrix, which includes 'G1' and 'G2' (VanRaden 2008). 
#' @return A n by n additive genomic relationship matrix. 
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' 
#' @example man/examples/computeG.R
#' 
#' @references \emph{VanRaden, P.M., 2008. Efficient methods to compute genomic predictions. Journal of dairy science, 91(11), pp.4414-4423.}
#' @export
computeG <- function(snpmatrix, maf, impute = 'mean', method) {
  if(anyNA(snpmatrix)) {
    for(j in 1 : ncol(snpmatrix)) {
      if(impute == 'mean') {
        snpmatrix[, j] <- ifelse(is.na(snpmatrix[, j]), mean(snpmatrix[, j], na.rm = TRUE), snpmatrix[, j])
      } else if(impute == 'rbinom') {
        set.seed(007) # reproducible imputation
        p_temp <- (colMeans(snpmatrix, na.rm = T) / 2)
        snpmatrix[, j] <- ifelse(is.na(snpmatrix[, j]), rbinom(1, 2, p_temp[j]), snpmatrix[, j])
      }
    }
  }
  p <- colMeans(snpmatrix) / 2
  maf_df <- pmin(p, 1-p)
  maf.index <- which(maf_df < maf)
  W <- snpmatrix[, -maf.index]
  if (method == 'G1') { 
    W_c <- scale(W, center = TRUE, scale = FALSE)
    G1 <- tcrossprod(W_c) / sum(2 * p * (1 - p))
    return(G1)
  } else if (method == 'G2'){
    W_cs <- scale(W, center = TRUE, scale = TRUE)
    G2 <- tcrossprod(W_cs) / ncol(W_cs)
    return(G2)
  }
}