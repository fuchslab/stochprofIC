#' Run the screening phase
#'
#' Given a tidy data frame with variables identifying the sample, gene, type, covariate(s) and 
#' transcripts per milion (TPM), run a pipeline based on a generalized additive model for 
#' location, shape and scale (GAMLSS) to detect/parameterize transcriptional heterogeneity 
#' of 10-cell data, normalizing for a set of covariates and adjusting for technical variation,
#' in order to obtain a data frame with the corrected data. Two consecutive GAMLSS are performed:
#' the first based on gene expression (expr) that is non-log, trimmed and zeroes converted to NAs, 
#' and the second based on log expr+1 without trimming.
#'
#' @encoding UTF-8
#' @param d Tidy data frame containing information about small-pool measurements. Each row corresponds to one measurement. Columns contain measured gene expression (\code{expr}), gene names (\code{gene}), measurement type (\code{type}) and covariate information to be adjusted for.
#' @param expr Name (character string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million (TPM). This column is expected to contain double values.
#' @param gene Name (character string) of the column in \code{d} that contains the names of the genes. This column is expected to contain factor values.
#' @param type Name (character string) of the column in \code{d} that contains the information whether the observation is from stochastic profiling (small-pool) samples or pool-and-split controls. This column is expected to contain factor values with exactly two levels: the pool-and-split controls and the small-pool samples.
#' @param covariates Names (vector of character strings) of the columns in \code{d} that contain further information to be adjusted for, aside of \code{type}. These columns are expected to contain factor values such as sex, biopsy, time, etc.
#' @param log.1st Logical value indicating whether a log transformation should be performed in the first GAMLSS. Defaults to FALSE.
#' @param log.2nd Logical value indicating whether a log transformation should be performed in the second GAMLSS. Defaults to TRUE.
#' @param add1.1st Logical value indicating whether a unit should be added to expr (expr+1) in the first GAMLSS. Defaults to FALSE.
#' @param add1.2nd Logical value indicating whether a unit should be added to expr (expr+1) in the second GAMLSS. Defaults to TRUE.
#' @param trim.1st Numerical value indicating the number of values to be trimmed, depending on each status of technical variation and covariates. When NA, no trimming is performed. Defaults to 3 in the first GAMLSS.
#' @param trim.2nd Numerical value indicating the number of values to be trimmed, depending on each status of technical variation and covariates. When NA, no trimming is performed. Defaults to NA in the second GAMLSS.
#' @param zeros.1st Logical value indicating whether to change the value of zero in expr, in the first GAMLSS.
#' @param zeros.2nd Logical value indicating whether to change the value of zero in expr, in the second GAMLSS.
#' @param zero.value.1st Numerical value with the substitute value for zero in expr, in the first GAMLSS.
#' @param zero.value.2nd Numerical value with the substitute value for zero in expr, in the second GAMLSS.
#' @param re Character string with the name of the variable that carries random effects. Not yet implemented.
#' @param re.placement Character string with either "location", "shape" or "both" as to whether the random effects in \code{re} are estimated for each of the components, or for both of them. Not yet implemented.
#' @export
#' @return A tibble with the screening values
#' @examples
#'\dontrun{
#'data(glioma)
#'stochprofIC.screening(glioma)
#'}
stochprofIC.screening <- function(...) {
  stochprofIC.pipeline(..., return = "screening")
}


