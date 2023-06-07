#' Run a gamlss() call in a specific way suitable for pipeline()
#'
#' Internal function to facilitate the call to gamlss() using the necessary
#' parts for pipeline()
#'
#' @param d Tidy data frame that has already been processed by pipeline() and contains a single gene.
#' @param return Character vector describing the type of output wanted: either "simple" to get the parameter estimates or "complete" to also include the corrected values
#' @param log See \code{\link{pipeline}}
#' @param trim See \code{\link{pipeline}}
#' @param zeros See \code{\link{pipeline}}
#' @param zero.value See \code{\link{pipeline}}
#' @param add.1 See \code{\link{pipeline}}
#' @param re See \code{\link{pipeline}}. Not yet implemented.
#' @param re.placement See \code{\link{pipeline}}. Not yet implemented.
#' @return Depending of the parameter "return", the output can be a tibble of parameter estimates ("simple") or also another tibble with the corrected values ("full").
run.gamlss <- function(d = d,
                       return = "simple",
                       log = NULL, trim = NULL, 
                       zeros = NULL, zero.value = NULL, add.1 = NULL,
                       re = NA, re.placement = "both") {

  # Deal with trimming the data
  # by removing the specific number of observations
  # in both technical variation 'tv' groups
  if (!is.na(trim)) {
    d <- d %>%
      dplyr::group_by(dplyr::across(-c(id, gene, TPM))) %>% 
      dplyr::arrange(TPM) %>%
      dplyr::slice(-c(1:trim)) %>%
      dplyr::arrange(dplyr::desc(TPM)) %>% 
      dplyr::slice(-c(1:trim)) %>%
      dplyr::ungroup()
  }
  # Deal with zeros
  if (zeros == TRUE) {
    d <- d %>%
      dplyr::mutate(TPM = ifelse(TPM == 0, zero.value, TPM))
  }
  # In any case, keep missing values out
  d <- dplyr::filter(d, !is.na(TPM))
  # Add 1 to TPM (tpm+1)
  if (add.1 == TRUE) {
    d <- d %>%
      dplyr::mutate(TPM = TPM + 1)
  }
  # Log or not the gene expression
  if (log == TRUE) {
    d <- d %>%
      dplyr::mutate(TPM = log(TPM))
  }

  # Define the equation
  model.equation.mu <- formula(paste0("TPM ~ ", paste(names(d)[!(names(d) %in% c("TPM", "id", "gene"))], 
                                                      collapse = " + ")))
  model.equation.sigma <- formula(paste0("~", paste(names(d)[!(names(d) %in% c("TPM", "id", "gene"))], 
                                                    collapse = " + ")))

  m <- invisible(try(gamlss::gamlss(model.equation.mu, sigma.formula = model.equation.sigma, 
                            data = d, family = NO, 
                            control = gamlss::gamlss.control(trace=FALSE))))

  if(return == "simple") {
    if(inherits(m, "try-error")){
      return(message("ERROR"))
    } else{
      sm <- invisible(try(summary(m)))
      if(inherits(sm, "try-error")){
        return(message("ERROR"))
      } else {
        out <- sm %>%
          tibble::tibble(rownames = "Variable") %>%
          dplyr::slice(((dplyr::n()/2)+1):dplyr::n()) %>%
          dplyr::select(Variable, Estimate, Pval = 5) %>%
          dplyr::filter(Variable == "tv")
        if(dim(out)[1] > 0){
          return(out)
        } else{
          return(tibble::tibble(Variable = "tv", Estimate = NA, Pval = NA))
        }
      }
    }
  } else{
    sm <- invisible(try(summary(m)))
    if(inherits(sm, "try-error")){
      return(message("ERROR"))
    } else {
      out <- sm %>%
        tibble::tibble(rownames = "Variable") %>%
        dplyr::mutate(Part = gl(2,(dplyr::n()/2), labels = c("Mu", "Sigma"))) %>%
        dplyr::select(Variable, Part, Estimate)
      mu.vector <- out %>%
        dplyr::filter(Part == "Mu") %>%
        dplyr::select(Variable, Estimate) %>%
        tibble::column_to_rownames("Variable") %>%
        as.matrix()
      sigma.vector <- out %>%
        dplyr::filter(Part == "Sigma") %>%
        dplyr::select(Variable, Estimate) %>%
        dplyr::mutate(Estimate = ifelse(Variable == "tv", 0, Estimate)) %>%
        tibble::column_to_rownames("Variable") %>%
        as.matrix()
      # d contains id, gene, TPM, and then the tv and the remaining covariates
      # Therefore the matrix multiplication must reflect that
      mu.part <- mu.vector[1] + (as.matrix(d[,-c(1,2,3)])%*% mu.vector[-1])
      sigma.part <- sigma.vector[1] + (as.matrix(d[,-c(1,2,3)]) %*% sigma.vector[-1])
      correction <- exp(((d$TPM - mu.part) / exp(sigma.part)))[,1]
      return(list(coefficients = out, TPMs = tibble::tibble(`Non-corrected` = d$TPM[d$tv == 1],
                                                    Corrected = correction[d$tv == 1])))
    }
  }
}


#' Extract Sidak Corrected p-values as a simple vector
#'
#' @param x a collection of p-values to be corrected using the Sidak correction
#' @return A vector of Sidak-corrected p-values.
extract.sidak <- function(x){
  return(mutoss::sidak(x)$adjPValues)
}


#' Extract Bonferroni Corrected p-values as a simple vector
#'
#' @param x a collection of p-values to be corrected using the Bonferroni correction
#' @return A vector of Bonferroni-corrected p-values.
extract.bonferroni <- function(x){
  return(mutoss::bonferroni(x)$adjPValues)
}

