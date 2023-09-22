#' Extract Sidak Corrected p-values as a simple vector
#'
#' @param x a collection of p-values to be corrected using the Sidak correction
#' @return A vector of Sidak-corrected p-values.
extract.sidak <- function(x){
  return(mutoss::sidak(x)$adjPValues)
}
