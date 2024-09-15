#' Run the whole pipeline
#'
#' Given a tidy data frame with variables identifying the sample, gene, type, covariate(s) and 
#' transcripts per milion (TPM), run a pipeline based on a generalized additive model for 
#' location, shape and scale (GAMLSS) to detect and parameterize transcriptional heterogeneity 
#' of 10-cell data, normalizing for a set of covariates and adjusting for technical variation,
#' to obtain the corrected data. Two consecutive GAMLSS are performed:
#' the first based on gene expression (expr) that is non-log, trimmed and zeroes converted to NAs, 
#' and the second based on log expr+1 without trimming.
#'
#' @encoding UTF-8
#' @param d Tidy data frame containing information about small-pool measurements. Each row corresponds to one measurement. Columns contain measured gene expression (\code{expr}), gene names (\code{gene}), measurement type (\code{type}) and covariate information to be adjusted for.
#' @param expr Name (character string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million (TPM). This column is expected to contain double values.
#' @param gene Name (character string) of the column in \code{d} that contains the names of the genes. This column is expected to contain factor values.
#' @param type Name (character string) of the column in \code{d} that contains the information whether the observation is a 10-cell biological sample or a 10-cell pool-and-split control. This column is expected to contain factor values with exactly two levels: 1) sample, 2) control.
#' @param covariates Names (vector of character strings) of the columns in \code{d} that contain further information to be adjusted for, aside from \code{type}. These columns are expected to contain factor values such as sex, biopsy, time, etc.
#' @param log.1st Logical value indicating whether a log transformation should be performed in the first GAMLSS. Defaults to FALSE.
#' @param log.2nd Logical value indicating whether a log transformation should be performed in the second GAMLSS. Defaults to TRUE.
#' @param add1.1st Logical value indicating whether a unit should be added to expr (expr+1) in the first GAMLSS. Defaults to FALSE.
#' @param add1.2nd Logical value indicating whether a unit should be added to expr (expr+1) in the second GAMLSS. Defaults to TRUE.
#' @param trim.1st Numerical value indicating the number of values to be trimmed, depending on each status of technical variation and covariates. When NA, no trimming is performed. Defaults to 3 in the first GAMLSS.
#' @param trim.2nd Numerical value indicating the number of values to be trimmed, depending on each status of technical variation and covariates. When NA, no trimming is performed. Defaults to NA in the second GAMLSS.
#' @param zeros.1st Logical value indicating whether to change the value of zero in expr, in the first GAMLSS.
#' @param zeros.2nd Logical value indicating whether to change the value of zero in expr, in the second GAMLSS.
#' @param zero.value.1st Numerical value with the substitute value for zero in expr, in the first GAMLSS.
#' @param zero.value.2nd Numerical value with the substitute value for zero in expr, in the second GAMLSS.
#' @param re Character string with the name of the variable that carries random effects. Not yet implemented.
#' @param re.placement Character string with either "location", "shape" or "both" as to whether the random effects in \code{re} are estimated for each of the components, or for both of them. Not yet implemented.
#' @param ... Other options, mainly the details about what to return. Used internally for \code{\link{stochprofIC.screening()}} and \code{\link{stochprofIC.testing}}.
#' @export
#' @return A list with two elements: "Coefficients" with a tibble containing gene, Variable, Part (either location "Mu" or shape "Sigma"), and Estimate; and "Corrected output" with a tibble containing gene and the non-corrected and corrected values.
#' @examples
#'\dontrun{
#'data(glioma)
#'stochprofIC.pipeline(glioma)
#'}
stochprofIC.pipeline <- function(d,
                            expr,
                            gene,
                            type,
                            covariates = NULL,
                            re = NA, re.placement = "both",
                            log.1st = FALSE, log.2nd = TRUE,
                            add1.1st = FALSE, add1.2nd = TRUE,
                            trim.1st = 3, trim.2nd = NA,
                            zeros.1st = TRUE, zeros.2nd = FALSE,
                            zero.value.1st = NA, zero.value.2nd = 0,
                            return = "Correcting&Inferring") {

  D <- arrange.input(d = d, expr = expr, gene = gene, type = type, covariates = covariates)
  N.total.genes <- length(unique(D$gene))

  # Keep a record of genes that give errors
  error.genes <- NULL
    
  ###########################################       
  ## Run the chosen steps of the pipeline: 
  ##   i) screening
  ##  ii) testing
  ## iii) correcting
  ##  iv) inferring
  # I suggest to combine iii) and iv) into a single "correcting&inferring" output
  ###########################################       
  
  ###########################################       
  ## Step i: screening
  ## 
  ## Estimate a GAMLSS;
  ## keep only those genes where samples have larger technical variation than controls;
  ## and where the p-value of the sigma coefficient for the technical variation (sigma.type) 
  ## is <= 0.05 after Sidak correction
  ## 
  ###########################################
      
  #-----------------------------------------------------
  # First call to 'stochprofIC.gamlss()'
  # do it gene by gene in a loop, and put it all back together into the object called 'out.gamlss1', 
  # that will contain a structured dataframe with the minimal required results 
  # (gene, sigma coefficient for technical variation, p-value for technical variation)
  #
  out.gamlss1 <- NULL
  list.genes <- levels(D$gene)
  names.covariates <- names(D)[names(D) %in% c("id", "gene", "expr", "type")]
  for (g in 1:length(list.genes)) {
    gene.now <- list.genes[g]
    d.gene <- D %>%
      filter(gene == gene.now)
    go <- NULL
    go <- stochprofIC.gamlss(d.gene,
                             expr = "expr", type = "type", gene = "gene", covariates = names.covariates,
                             re = re, re.placement = re.placement,
                             log = log.1st,
                             add.1 = add1.1st,
                             trim = trim.1st, 
                             zeros = zeros.1st, 
                             zero.value = zero.value.1st,
                             input.cleaned = TRUE)
    if (is.null(go)) {
      error.genes <- c(error.genes, gene.now)
    } else {
      out.gamlss1 <- bind_rows(out.gamlss1, 
                     bind_cols(gene = gene.now, go))
    }
  }

  
  
  #-----------------------------------------------------
  # Second part: Run a Sidak Correction of sigma-type p-values and keep only
  # sigma.type > 0 & Pval.sidak <= 0.05
  
  out.sidak <- out.gamlss1 %>%
    dplyr::mutate(Pval.sidak = extract.sidak(Pval))
  
  out.screening1 <- out.gamlss1 %>%
    dplyr::mutate(Pval.sidak = extract.sidak(Pval)) %>%
    dplyr::filter(Estimate >= 0) %>%
    dplyr::filter(Pval.sidak <= 0.05)
  
  if(return == "screening"){
    # Proportion of genes kept after sidak correction and removing negative sigma estimates
    message(paste("Proportion of genes kept after sidak correction\nand removing negative sigma estimates:", (dim(out.screening1)[1] / N.total.genes)*100, "%"))
    return(list(Screening = out.screening1, 
                `Genes with error` = unique(error.genes)))
  }

  
  ###########################################       
  ## Step ii: testing
  ## 
  ## Estimate a GAMLSS to identify possible batch effects
  ## 
  ###########################################

  
  #-----------------------------------------------------
  # Third part: Run a GAMLSS with log of expr+1.
  # Return a single row for each gene with the Sigma coefficient and sigma p-value
  # for the variable that captures the technical/biological variation (type) as well as
  # the rest of mu and sigma coefficients to correct for technical variation
  # It requires 
  
  # Run several genes
  d.post.sidak <- D %>%
    dplyr::filter(gene %in% out.screening1$gene) %>%
    dplyr::mutate(gene = factor(gene)) %>% 
    droplevels()
  out.gamlss2 <- NULL
  for(g in 1:length(levels(d.post.sidak$gene))){
    gene.now <- levels(d.post.sidak$gene)[g]
    go <- NULL
    go <- stochprofIC.gamlss(dplyr::filter(d.post.sidak, 
                       gene == levels(d.post.sidak$gene)[g]),
                       expr = "expr", type = "type", gene = "gene",
                       re = re, re.placement = re.placement,
                       log = log.2nd,
                       add.1 = add1.2nd,
                       trim = trim.2nd, 
                       zeros = zeros.2nd, 
                       zero.value = zero.value.2nd,
                       return = "simple",
                       input.cleaned = TRUE)
    if (is.null(go)) {
      error.genes <- c(error.genes, gene.now)
    } else {
      out.gamlss2 <- dplyr::bind_rows(out.gamlss2, 
            dplyr::bind_cols(gene = levels(d.post.sidak$gene)[g], go))
    }
  }
  
  #-----------------------------------------------------
  # Fourth part: Run a FDR Correction of sigma-type p-values and keep only
  # sigma.type > 0 & Pval.fdr <= 0.05
  
  # Manual test for extracting bonferroni p-values (FDR)
  out.fdr <- out.gamlss2 %>%
    dplyr::mutate(Pval.dfr = stats::p.adjust(Pval, method="fdr")) 
  
  out.screening2 <- out.gamlss2 %>%
    dplyr::mutate(Pval.fdr = stats::p.adjust(Pval, method="fdr")) %>%
    dplyr::filter(Estimate >= 0) %>%
    dplyr::filter(Pval.fdr <= 0.05)
  
  
  
  if(return == "testing"){
    # Proportion of genes kept after FDR correction and removing negative sigma estimates
    message(paste("Proportion of genes kept after FDR correction\nand removing negative sigma estimates:", (dim(out.screening2)[1] / N.total.genes)*100, "%"))
    return(list(Testing = out.screening2, 
                `Genes with error` = unique(error.genes)))
  }
  

  ###########################################       
  ## Step iii: correcting
  ## 
  ## Use the estimated GAMLSS coefficients from the previous testing step
  ## for batch correction
  ## (this only makes sense if testing was done)
  ## 
  ###########################################

  
  
  #----------------------------------------------------
  # Fifth part: Run again the GAMLSS with log of expr+1 only for selected genes.
  # Return a full corrected values for technical variation only in Samples
  
  d.post.fdr <- D %>%
    dplyr::filter(gene %in% out.screening2$gene) %>%
    droplevels()

  
  out.corrected <- out.coefficients <- NULL
  for(g in 1:length(levels(d.post.fdr$gene))){
    d.now <- d.post.fdr %>%
      dplyr::filter(gene == levels(d.post.fdr$gene)[g])
    out.full <- NULL
    out.full <- stochprofIC.gamlss(d.now,
                                   expr = "expr", type = "type", gene = "gene",
                           return = "full", 
                           re = re, re.placement = re.placement,
                           log = log.2nd,
                           add.1 = add1.2nd,
                           trim = trim.2nd, 
                           zeros = zeros.2nd, 
                           zero.value = zero.value.2nd,
                           input.cleaned = TRUE)

    if (is.null(out.full)) {
      error.genes <- c(error.genes, levels(d.post.fdr$gene)[g])
    } else {
      out.corrected <- dplyr::bind_rows(out.corrected, 
                                 tibble::tibble(gene = levels(d.post.fdr$gene)[g],
                                        out.full$exprs))   # Christiane: maybe change the name "exprs" in the gamlss function -- a bit confusing
      out.coefficients <- dplyr::bind_rows(out.coefficients, 
                                    tibble::tibble(gene = levels(d.post.fdr$gene)[g],
                                    out.full$coefficients))
    }
  }
  
  out.corrected <- out.corrected %>%
    dplyr::mutate(gene = factor(gene))
  out.coefficients <- out.coefficients %>%
    dplyr::mutate(gene = factor(gene),
           Variable = factor(Variable))
  
  # Columns have to be renamed again!  
  
  # Screening information: 
  # Proportion of genes kept after sidak correction and removing negative sigma estimates
  message(paste("Proportion of genes kept after sidak correction\nand removing negative sigma estimates:", (dim(out.screening1)[1] / N.total.genes)*100, "%"))
  # Testing information: 
  # Proportion of genes kept after FDR correction and removing negative sigma estimates
  message(paste("Proportion of genes kept after FDR correction\nand removing negative sigma estimates:", (dim(out.screening2)[1] / N.total.genes)*100, "%"))
  
  return(list(Coefficients = out.coefficients,
              `Corrected output` = out.corrected,
              `Genes with error` = unique(error.genes)))
}


