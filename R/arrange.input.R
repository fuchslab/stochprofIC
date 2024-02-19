#' Arrange the raw input and tidy it for both the pipeline and the gamlss standalone functions
#'
#' Internal function to arrange the raw input and arrange it into known variable names, ready for the analysis.
#'
#' @param d Tidy data frame containing information about small-pool measurements. Each row corresponds to one measurement. Columns contain measured gene expression (\code{expr}), gene names (\code{gene}), measurement type (\code{type}) and covariate information to be adjusted for.
#' @param expr Name (character string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million (TPM). This column is expected to contain double values.
#' @param gene Name (character string) of the column in \code{d} that contains the names of the genes. This column is expected to contain factor values.
#' @param type Name (character string) of the column in \code{d} that contains the information whether the observation is from stochastic profiling (small-pool) samples or pool-and-split controls. This column is expected to contain factor values with exactly two levels: the pool-and-split controls and the small-pool samples.
#' @param covariates Names (vector of character strings) of the columns in \code{d} that contain further information to be adjusted for, aside of \code{type}. These columns are expected to contain factor values such as sex, biopsy, time, etc.
#' @return A tidy dataset with the data arranged, cleaned and ready for the analysis.
arrange.input <- function(d, expr, gene, type, covariates = NULL) {

  #############################################################
  ## Rearranging input matrix and checking of correct format ##
  #############################################################
  
  # Add an additional column "id", assigning an internal ID to each row
  d <- dplyr::mutate(d, id = 1:dplyr::n())
  
  # Ensure the column containing the gene expression values is called "expr"
  d <- dplyr::rename(d, expr = !!expr)
  
  # Ensure that "expr" contains numerical value
  if (!is.double(d$expr)) {
    stop("'expr' must contain values of type double.")
  }
  
  # Ensure the column containing the gene names is called "gene"
  d <- dplyr::rename(d, gene = !!gene)
  
  # Ensure the column containing the measurement type is called "type"
  d <- dplyr::rename(d, type = !!type)


# Christiane: Code checking to be continues from here

  # Parser technical variation, ensure that it is a factor with only 2 levels
  # and that its name is 'type'
  if (!is.factor(d$type)) {
    stop("'type' must be the name of an existing factor, 
         whose first level refers to the pool-and-split controls.")
  }
  if (length(levels(d$type)) != 2) {
    stop("The technical variation must be a factor with strictly 2 levels, 
         the first one referring to the pool-and-split controls, 
         and the second one to the small-pool samples.")
  }

  # Ensure that for technical variation (type) zero and one are created
  # and that the resulting value is an integer (with 1L)
  d <- d %>% dplyr::mutate(type = as.integer(type) - 1L)

  # Parse the covariates by first selecting them 
  # and separating them in a different object named 'd.cv'
  d.cv <- dplyr::select(d, id, tidyselect::all_of(covariates))
  # Ensure that all of them are factors
  class.covariates <- unlist(lapply(dplyr::select(d.cv, -id), class))
  if (length(class.covariates) != length(which(class.covariates %in% "factor"))) {
    stop("All variables in 'covariates' must be factors.")
  }
  # Create a model matrix with all covariates, and exclude the intercept
  d.mm <- stats::model.matrix(~ ., data = dplyr::select(d.cv, -id)) %>%
    data.frame() %>%
    select(-1)
  # Keep only the relevant variables in a tidy format
  # From now on, 'D' is the result of arranging and checking the input
  # provided by the user. It contains the variables in a specific order: id, gene, expr, type, covariates.
  D <- dplyr::select(d, id, gene, expr, type) %>%
    dplyr::left_join(dplyr::mutate(d.mm, id = 1:dplyr::n()), by = "id")
  # 
  return(D)

}

