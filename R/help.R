#' Tools facilitating a pipeline based on a generalized additive model to detect/parameterize transcriptional heterogeneity of 10-cell data.
#'
#' Development of a pipeline to detect/parameterize transcriptional heterogeneity of 10-cell data, normalizing for a set of covariates and adjusting for technical variation, in order to obtain a the corrected data.
#'
#' @encoding UTF-8
#' @references \url{http://mercegari.net/}.
#' @importFrom dplyr %>% mutate filter select full_join tibble rename group_by arrange across bind_cols bind_rows ungroup n desc 
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist NO
#' @importFrom mutoss sidak bonferroni
#' @importFrom stats formula model.matrix p.adjust
#' @importFrom tidyselect all_of
#' @docType package
#' @name pipeline
NULL
