#' Run a gamlss() call in a specific way suitable for stochprofIC.pipeline
#'
#' Internal function to facilitate the call to gamlss() using the necessary
#' parts for stochprofIC.pipeline
#'
#' @param d Tidy data frame that has already been processed by stochprofIC.pipeline and contains a single gene.
#' @param return Character vector describing the type of output wanted: either "simple" to get the parameter estimates or "complete" to also include the corrected values
#' @param log See \code{\link{stochprofIC.pipeline}}
#' @param trim See \code{\link{stochprofIC.pipeline}}
#' @param zeros See \code{\link{stochprofIC.pipeline}}
#' @param zero.value See \code{\link{stochprofIC.pipeline}}
#' @param add.1 See \code{\link{stochprofIC.pipeline}}
#' @param re See \code{\link{stochprofIC.pipeline}}. Not yet implemented.
#' @param re.placement See \code{\link{stochprofIC.pipeline}}. Not yet implemented.
#' @return Depending of the parameter "return", the output can be a tibble of parameter estimates ("simple") or also another tibble with the corrected values ("full").
stochprofIC.gamlss <- function(d,
                       return = "simple",
                       log = NULL, trim = NULL, 
                       zeros = NULL, zero.value = NULL, add.1 = NULL,
                       re = NA, re.placement = "both") {
  
# The following lines assume that the column names are exactly as defined in stochprofIC.pipeline. We should not do that.
  ## I don't understand why not... I guess there are 2 issues to discuss here:
  ## First, do we really want this function to work independently of the stochprofIC.pipeline? If not, the rest of the discussion is over.
  ## I personally think that there is no need to have this function exported as a general function in the package. 
  ## However (and second), if we want to have this stochprofIC.gamlss as a general function in the package, then certainly
  ## we need to ensure that the input object contains the names that we want. I can improve the documentation (line 6).
  
  # Deal with trimming the data
  # by removing the specific number of observations
  # in both technical variation 'type' groups
  if (!is.na(trim)) {
    d <- d %>%
      dplyr::group_by(dplyr::across(-c(id, gene, expr))) %>% 
      dplyr::arrange(expr) %>%
      dplyr::slice(-c(1:trim)) %>%
      dplyr::arrange(dplyr::desc(expr)) %>% 
      dplyr::slice(-c(1:trim)) %>%
      dplyr::ungroup()
  }
  # Deal with zeros
  if (zeros == TRUE) {
    d <- d %>%
      dplyr::mutate(expr = ifelse(expr == 0, zero.value, expr))
  }
  # In any case, keep missing values out
  d <- dplyr::filter(d, !is.na(expr))
  # Add 1 to expr (expr+1)
  if (add.1 == TRUE) {
    d <- d %>%
      dplyr::mutate(expr = expr + 1)
  }
  # Log or not the gene expression
  if (log == TRUE) {
    d <- d %>%
      dplyr::mutate(expr = log(expr))
  }
  
  # Build the model equations for the GAMLSS
  
  # Equation for the location parameter:
  # on the left expr,
  # on the right all column names except expr, id and gene (each as a linear term)
  
  model.equation.mu <- formula(paste0("expr ~ ", paste(names(d)[!(names(d) %in% c("expr", "id", "gene"))], 
                                                      collapse = " + ")))
  
  # Equation for the scale parameter:
  # on the right all column names except expr, id and gene (each as a linear term)
  model.equation.sigma <- formula(paste0("~", paste(names(d)[!(names(d) %in% c("expr", "id", "gene"))], 
                                                    collapse = " + ")))
  
  # Estimation of the GAMLSS
  m <- invisible(try(gamlss::gamlss(formula = model.equation.mu, 
                                    sigma.formula = model.equation.sigma, 
                                    data = d, family = gamlss.dist::NO, 
                                    control = gamlss::gamlss.control(trace=FALSE))))
  
  
## I (CF) don't know what the folloing lines are doing
  
## Why try and error here?  
  ## Try and error allows us to continue with the function when gamlss has failed
  ## or not converged or not produced all the necessary estimates.
  
  if(return == "simple") {
    if(inherits(m, "try-error")){
      return(message("ERROR"))
    } else{
      sm <- invisible(try(summary(m)))
      if(inherits(sm, "try-error")){
        return(message("ERROR"))
      } else {
        out <- sm %>%
          as.data.frame() %>% 
          tibble::rownames_to_column("Variable") %>%
          dplyr::select(Variable, Estimate, Pval = 5) %>%
          dplyr::filter(Variable == "type.1")
        if(dim(out)[1] > 0){
          return(out)
        } else{
          return(tibble::tibble(Variable = "type", Estimate = NA, Pval = NA))
        }
      }
    }
  } else{
    sm <- invisible(try(summary(m)))
    if(inherits(sm, "try-error")){
      return(message("ERROR"))
    } else {
      out <- sm %>%
        as.data.frame() %>% 
        tibble::rownames_to_column("Variable") %>%
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
        dplyr::mutate(Estimate = ifelse(Variable == "type", 0, Estimate)) %>%
        tibble::column_to_rownames("Variable") %>%
        as.matrix()
  
# The following refers to fixed column numbers 1,2,3. This is too error-prone because the (ir)relevant columns might be anywhere. --> Change
      ## In principle, we ensure the structure of d in the stochprofIC.pipeline.R, line 117.
      ## Therefore, the order will be always 1=id, 2=gene and 3=expr.
# Correction is already part of this function (--> function should have a different name and description).
      ## Again, this depends on what we want about stochprofIC.gamlss, if it is a function itself or not.
      ## But if it has to be out, this is a big design change, because from gamlss I already 
      ## have both original data and gamlss estimates, while in the stochprofIC.pipeline, I don't 
      ## have both in the same way. I could try...
# If correction happens here -- what happens to those values that were trimmed?
      ## In the stochprofIC.pipeline function we already define the trimming 
      ## for the 1st gamlss and for the 2nd gamlss: trim.1st = 3, trim.2nd = NA. 
      ## Therefore, in the 2nd gamlss, in which we should not have trimmed,
      ## we obtain the corrected values for all the input data.
# Is log-transformation taken into account?
      ## The log-transformation is only taken into account in the 2nd gamlss.
      ## By default, the stchorpfoIC.pipeline has: log.1st = FALSE, log.2nd = TRUE
      
      # d contains id, gene, expr, and then the type and the remaining covariates
      # Therefore the matrix multiplication must reflect that
      mu.part <- mu.vector[1] + (as.matrix(d[,-c(1,2,3)])%*% mu.vector[-1])
      sigma.part <- sigma.vector[1] + (as.matrix(d[,-c(1,2,3)]) %*% sigma.vector[-1])
      correction <- exp(((d$expr - mu.part) / exp(sigma.part)))[,1]
      
      return(list(coefficients = out, exprs = tibble::tibble(`Non-corrected` = d$expr[d$type == 1],
                                                            Corrected = correction[d$type == 1])))
    }
  }
}
