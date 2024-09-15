#' Run the testing phase
#'
#' Given a tidy data frame with variables identifying the sample, gene, type, covariate(s) and 
#' transcripts per milion (TPM), test a generalized additive model for 
#' location, shape and scale (GAMLSS) to detect and parameterize transcriptional heterogeneity 
#' of 10-cell data, normalizing for a set of covariates and adjusting for technical variation,
#' to obtain a data frame of tested genes. Two consecutive GAMLSS are performed:
#' the first based on gene expression (expr) that is non-log, trimmed and zeroes converted to NAs (the same as \code{stochprofIC.screening()}), 
#' and the second based on log expr+1 without trimming. With this testing phase, a GAMLSS estimation is ready to identify possible batch effects.
#'
#' @encoding UTF-8
#' @param d Tidy data frame containing information about 10-cell measurements. Each row corresponds to one measurement. Columns contain measured gene expression (\code{expr}), gene names (\code{gene}), measurement type (\code{type}) and covariate information to be adjusted for.
#' @param expr Name (character string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million (TPM). This column is expected to contain double values.
#' @param gene Name (character string) of the column in \code{d} that contains the names of the genes. This column is expected to contain factor values.
#' @param type Name (character string) of the column in \code{d} that contains the information whether the observation is from a 10-cell biological sample or a 10-cell pool-and-split control. This column is expected to contain factor values with exactly two levels: 1) sample, 2) control.
#' @param covariates Names (vector of character strings) of the columns in \code{d} that contain further information to be adjusted for, aside from \code{type}. These columns are expected to contain factor values such as sex, biopsy, time, etc.
#' @param log Logical value indicating whether a log transformation should be performed in the testing GAMLSS. Defaults to TRUE.
#' @param add1 Logical value indicating whether a unit should be added to expr (expr+1) in the testing GAMLSS. Defaults to TRUE.
#' @param trim Numerical value indicating the number of values to be trimmed, depending on each status of technical variation and covariates. When NA, no trimming is performed. Defaults to NA in the testing GAMLSS.
#' @param zeros Logical value indicating whether to change the value of zero in expr, in the testing GAMLSS.
#' @param zero.value Numerical value with the substitute value for zero in expr, in the testing GAMLSS.
#' @param re Character string with the name of the variable that carries random effects. Not yet implemented.
#' @param re.placement Character string with either "location", "shape" or "both" as to whether the random effects in \code{re} are estimated for each of the components, or for both of them. Not yet implemented.
#' @export
#' @return A tibble with the testing values
#' @examples
#'\dontrun{
#'data(glioma)
#'stochprofIC.testing(glioma)
#'}
stochprofIC.testing <- function(d,
				expr,
				gene,
				type,
				covariates,
				log,
				add1,
				trim,
				zeros,
				zero.value,
				re = NA,
				re.placement = "both") {
  pipeline.args <- list(d = d,
		       expr = expr,
		       gene = gene,
		       type = type,
		       covariates = covariates,
		       log.2nd = log,
		       add1.2nd = add1,
		       trim.2nd = trim,
		       zeros.2nd = zeros,
		       zero.value.2nd = zero.value,
		       re = re,
		       re.placement = re.placement,
		       return = "testing")
  do.call(stochprofIC.pipeline, pipeline.args)
}

