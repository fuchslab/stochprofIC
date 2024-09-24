#' Run a customized gamlss() call for stochprofIC.pipeline
#'
#' Internal function to implement gamlss() using the necessary
#' parts for stochprofIC.pipeline
#'
#' @param d Tidy data frame for a single gene that has already been processed. If the data have not previously been arranged, the argument \code{is.input.cleaned} should be set to FALSE so that the input is arranged before analysis.
#' @param expr Name (character string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million (TPM). This column is expected to contain double values.
#' @param gene Name (character string) of the column in \code{d} that contains the names of the genes. This column is expected to contain factor values.
#' @param type Name (character string) of the column in \code{d} that contains the information whether the observation is from a 10-cell biological sample or a 10-cell pool-and-split control. This column is expected to contain factor values with exactly two levels: 1) sample, 2) control.
#' @param covariates Names (vector of character strings) of the columns in \code{d} that contain further information to be adjusted for, aside from \code{type}. These columns are expected to contain factor values such as sex, biopsy, time, etc.
#' @param return Character vector specifying the output: "simple" returns and "complete" appends the corrected values.
#' @param log Logical value indicating whether a log transformation should be performed.
#' @param trim Numerical value indicating the number of values to be trimmed, depending on each status of technical variation and covariates. When NA, no trimming is performed.
#' @param zeros Logical value indicating whether to change the value of zero in expr.
#' @param zero.value Numerical value with the substitute value for zero in expr.
#' @param add.1 Logical value indicating whether a unit should be added to expr (expr+1).
#' @param re Character string with the name of the variable that carries random effects. Not yet implemented.
#' @param re.placement Character string with either "location", "shape" or "both" as to whether the random effects in \code{re} are estimated for each of the components, or for both of them. Not yet implemented.
#' @param is.input.cleaned Logical Whether the input is already cleaned (when used within \code{\link{stochprofIC.pipeline}}). Defaults to FALSE, when the data is cleaned within.
#' @export
#' @return A tibble of parameter estimates (if "return" = "simple") or a tibble of parameter estimates and a second tibble of corrected values (if "return" = "full").
stochprofIC.gamlss <- function(d,
                       expr,
                       gene, 
                       type,
                       covariates = NULL,
                       return = "simple",
                       log = TRUE, trim = NA, 
                       zeros = FALSE, zero.value = 0, add.1 = TRUE,
                       re = NA, re.placement = "both",
                       is.input.cleaned = FALSE) {
  
  if (!is.input.cleaned) {
    d <- arrange.input(d = d, expr = expr, gene = gene, type = type, covariates = covariates)
  }
  gene.now <- unique(d$gene)

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
  capture.output({m <- try(suppressWarnings(gamlss::gamlss(formula = model.equation.mu, 
                          sigma.formula = model.equation.sigma, 
                          data = d, family = gamlss.dist::NO, 
                          control = gamlss::gamlss.control(trace=FALSE)))) }, file = nullfile())
  
  if(return == "simple") {
    if(inherits(m, "try-error")){
      return(message(paste0("  ERROR", " in gene ", gene.now)))
    } else{
      capture.output({ sm <- try(suppressWarnings(summary(m))) }, file = nullfile()) 
      if(inherits(sm, "try-error")){
        return(message(paste0("  ERROR", " in gene ", gene.now)))
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
    capture.output({ sm <- try(suppressWarnings(summary(m))) }, file = nullfile())
    if(inherits(sm, "try-error")){
      return(message(paste0("  ERROR", " in gene ", gene.now)))
    } else {
      keep.rownames <- rownames(sm)
      out <- sm %>%
        base::as.data.frame() %>% 
        dplyr::mutate(Variable = keep.rownames) %>% 
        dplyr::mutate(Part = gl(2,(dplyr::n()/2), labels = c("Mu", "Sigma"))) %>%
        dplyr::select(Variable, Part, Estimate) #%>%
        #dplyr::filter(!grepl('Intercept', Variable))
      rownames(out) <- NULL
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

      # Remove intercepts from the return
      out <- out %>% 
        dplyr::filter(!grepl('Intercept', Variable))
      
      return(list(coefficients = out, exprs = tibble::tibble(`Non-corrected` = d$expr[d$type == 1],
                                                            Corrected = correction[d$type == 1])))
    }
  }
}
