#' Run a gamlss() call in a specific way suitable for stochprofIC.pipeline
#'
#' Internal function to facilitate the call to gamlss() using the necessary
#' parts for stochprofIC.pipeline
#'
#' @param d Tidy data frame that has already been processed and that contains a single gene. If the data has not previously been arranged, the argument \code{input.cleaned} can be set to TRUE and then the input is arranged before the analysis.
#' @param expr Name (character string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million (TPM). This column is expected to contain double values.
#' @param gene Name (character string) of the column in \code{d} that contains the names of the genes. This column is expected to contain factor values.
#' @param type Name (character string) of the column in \code{d} that contains the information whether the observation is from stochastic profiling (small-pool) samples or pool-and-split controls. This column is expected to contain factor values with exactly two levels: the pool-and-split controls and the small-pool samples.
#' @param covariates Names (vector of character strings) of the columns in \code{d} that contain further information to be adjusted for, aside of \code{type}. These columns are expected to contain factor values such as sex, biopsy, time, etc.
#' @param return Character vector describing the type of output wanted: either "simple" to get the parameter estimates or "complete" to also include the corrected values
#' @param log See \code{\link{stochprofIC.pipeline}}
#' @param trim See \code{\link{stochprofIC.pipeline}}
#' @param zeros See \code{\link{stochprofIC.pipeline}}
#' @param zero.value See \code{\link{stochprofIC.pipeline}}
#' @param add.1 See \code{\link{stochprofIC.pipeline}}
#' @param re See \code{\link{stochprofIC.pipeline}}. Not yet implemented.
#' @param re.placement See \code{\link{stochprofIC.pipeline}}. Not yet implemented.
#' @param input.cleaned Logical Whether the input is already cleaned (when used within \code{\link{stochprofIC.pipeline}}). Defaults to FALSE, when the data is cleaned within.
#' @export
#' @return Depending on the parameter "return", the output can be a tibble of parameter estimates ("simple") or also another tibble with the corrected values ("full").
stochprofIC.gamlss <- function(d,
                       expr,
                       gene, 
                       type,
                       covariates = NULL,
                       return = "simple",
                       log = TRUE, trim = NA, 
                       zeros = FALSE, zero.value = 0, add.1 = TRUE,
                       re = NA, re.placement = "both",
                       input.cleaned = FALSE) {
  
  if (!input.cleaned) {
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
