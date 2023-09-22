#' Extract Bonferroni Corrected p-values as a simple vector
#'
#' @param x a collection of p-values to be corrected using the Bonferroni correction
#' @return A vector of Bonferroni-corrected p-values.
extract.bonferroni <- function(x){
  return(mutoss::bonferroni(x)$adjPValues)
}
