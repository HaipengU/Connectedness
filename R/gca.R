#' \strong{GCA}: Genetic connectedness analysis
#' 
#' @description 
#' An R package for genetic connectedness analysis across units using pedigree and genomic data.
#' 
#' @details 
#' This package implements three connectedness metrics which are functions of prediction error variance (PEV) matrix,  
#'  including prediction error variance of difference (PEVD), coefficient of determination (CD) and prediction error correlation (r), 
#'  coupled with three summary methods for each statistic. For example, the PEVD across units is summarized as 1) average PEV of all 
#'  pairwise differences between individuals across units; 2) average PEV within and across units; 3) using a contrast vector. 
#'  Analogous summary methods will be also applied to CD and r statistics. Three additional metrics approximating connectedness using 
#'  variance of estimates of unit effects (VE) are included, such as variance of estimates of units effects difference (VED), 
#'  coefficient of determination of VED (CDVED), and connectedness rating (CR). Within each metric, three different methods are names 
#'  according to the number of corrected factors (e.g., 0, 1 and 2) for fixed effects, such as VED0, VED1 and VED2. Similar corrected 
#'  function is also applied for metrics of CDVED and CR. 
#' 
#' @section Available functions in GCA package:
#' \itemize{
#'   \item computeA(): Computation of numerator relationship matrix. 
#'   \item computeG(): Computation of genomic relationship matrix.
#'   \item gcm(): Measures of genetic connectedness.
#'   \item varcomp(): Estimates of variance components using eigenvalues and eigenvectors. 
#' }
#'
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#' 
#' 
#'@keywords internal
"_PACKAGE"


