#' Cattle dataset 
#'
#' The cattle dataset is simulated with QMSim software (Sargolzaei and Schenkel, 2009).   
#' This dataset includes 530 individuals across five generations (from founder to generation 4), 
#' each with 10,000 single nucleotide polymorphisms spread over 29 autosomes. Single phenotype with heritability of 0.6 was simulated.
#' Two fixed covariates of sex and unit (management unit) are available.
#'
#' @docType data
#' @name GCcattle
#' @usage 
#' data(GCcattle)
#' head(cattle.pheno)
#' dim (cattle.W)
#' @format 
#' \itemize{
#'   \item cattle.pheno: contains phenotype, fixed covariates and pedigree information.
#'   \item catlle.W: contains marker information.
#' }
#'
#' @keywords datasets
#'
#' @references Sargolzaei, M., and F. S. Schenkel. 2009. Qmsim: a large-scale
#' genome simulator for livestock. Bioinformatics 25:680–681. doi:10.1093/bioinformatics/btp045
#'
NULL

