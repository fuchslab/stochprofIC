
#' Breast cancer: a dataset on  XX
#'
#' A dataset on X contining observations on X, ready for being used in the pipeline
#'
#' @name breast.cancer 
#' @docType data
#' @keywords datasets
#' @usage data(breast.cancer)
#' @format A data frame containing X suitable for being used in the pipeline.
#' @examples
#' data(breast.cancer)
#' str(breast.cancer)
#'\dontrun{
#'stochprofIC.pipeline(d = breast.cancer,
#'                expr = "TPM",
#'                gene = "gene", 
#'                type = "type",
#'                covariates = "biopsy",
#'                re = NA,   # not implemented yet
#'                log.1st = FALSE, 
#'                log.2nd = TRUE,
#'                add1.1st = FALSE, 
#'                add1.2nd = TRUE,
#'                trim.1st = 3, 
#'                trim.2nd = NA,
#'                zeros.1st = TRUE, zeros.2nd = FALSE,
#'                zero.value.1st = NA, zero.value.2nd = 0,
#'                re.placement = "both",
#'                return="screening")
#'}
NULL


#' Glioma: a dataset on  XX
#'
#' A dataset on X contining observations on X, ready for being used in the pipeline
#'
#' @name glioma
#' @docType data
#' @keywords datasets
#' @usage data(glioma)
#' @format A data frame containing X suitable for being used in the pipeline.
#' @examples
#' data(glioma)
#' str(glioma)
#'\dontrun{
#'stochprofIC.pipeline(d = glioma,
#'               expr = "TPM",
#'               type = "type",
#'               covariates = "sex",
#'               gene = "gene", 
#'               re = NA,
#'               log.1st = FALSE, 
#'               log.2nd = TRUE,
#'               add1.1st = FALSE, 
#'               add1.2nd = TRUE,
#'               trim.1st = 3, 
#'               trim.2nd = NA,
#'               zeros.1st = TRUE, zeros.2nd = FALSE,
#'               zero.value.1st = NA, zero.value.2nd = 0,
#'               re.placement = "both")
#'}
NULL

