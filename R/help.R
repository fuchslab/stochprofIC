#' Tools facilitating a pipeline based on a generalized additive model to detect and parameterize transcriptional heterogeneity of 10-cell data.
#'
#' A pipeline to detect and parameterize transcriptional heterogeneity of 10-cell data by accounting for a set of covariates and adjusting for technical variation to obtain corrected data.
#'
#' @encoding UTF-8
#' @references \url{http://mercegari.net/}.
#' @importFrom dplyr %>% mutate filter select full_join tibble rename group_by arrange across bind_cols bind_rows ungroup n desc 
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist NO
#' @importFrom mutoss sidak bonferroni
#' @importFrom stats formula model.matrix p.adjust
#' @importFrom utils capture.output
#' @importFrom tidyselect all_of
#' @name stochprofIC
"_PACKAGE"
#' @aliases stochprofIC
NULL
