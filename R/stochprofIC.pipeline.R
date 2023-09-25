#' Run the whole pipeline
#'
#' Given a tidy data frame with X, X, X, run a pipeline based on a generalized additive model
#' for location, shape and scale (GAMLSS) to detect/parameterize transcriptional heterogeneity 
#' of 10-cell data, normalizing for a set of covariates and adjusting for technical variation,
#' in order to obtain a data frame with the corrected data. Two consecutive GAMLSS are performed:
#' the first based on gene expression (expr) that is non-log, trimmed and zeroes converted to NAs, and the second based on 
#' log expr+1 without trimming.
#'
#' ### Merce, please check the description of the next lines for correctness.
#'
#' @encoding UTF-8
#' @param d Tidy data frame containing information about small-pool measurements. Each row corresponds to one measurement. Columns contain measured gene expression (\code{expr}), gene names (\code{gene}), measurement type (\code{type}) and covariate information to be adjusted for.
#' @param expr Name (string) of the column in \code{d} that contains the measurements of gene expression; e.g. as transcripts per million. This column is expected to contain double values.
#' @param gene Name (string) of the column in \code{d} that contains the names of the genes. This column is expected to contain string values.
#' @param type Name (string) of the column in \code{d} that contains the information whether the observation is from stochastic profiling (small-pool) samples or pool-and-split controls. This column is expected to contain factor values with exactly two levels.
#' @param covariates Names (vector of strings) of the columns in \code{d} that contain further information to be adjusted for, aside of \code{type}. These columns are expected to contain factor variables.
#' 
#' Needed: description of these two levels
#' 
#' ### stop checking here (I have to continue later)
#' 
#' 
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
#' @export
#' @return A list with two elements: "Coefficients" with a tibble containing gene, Variable, Part (either location "Mu" or shape "Sigma"), and Estimate; and "Corrected output" with a tibble containing gene and the non-corrected and corrected values.
#' @examples
#'\dontrun{
#'stochprofIC.pipeline()
#' XXXX we need to think a way to pass some example data
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
                    # the following functionality was discussed but is not yet properly implemented, not tested and not explained!
                    # the step "inferring" is missing completely
                            screening=T,testing=T,correcting=T,inferring=T) {

  #############################################################
  ## Rearranging input matrix and checking of correct format ##
  #############################################################
  
  # Add an additional column "id", assigning an internal ID to each row
  d <- dplyr::mutate(d, id = 1:dplyr::n())
  
  # Make sure the column containing the gene expression values is called "expr"
  d <- dplyr::rename(d, expr = !!expr)
  
    # Make sure the column containing the gene names is called "gene"
  d <- dplyr::rename(d, gene = !!gene)
  
  # Make sure the column containing the measurement type is called "type"
  d <- dplyr::rename(d, type = !!type)
  
  
  # Parser for expr: ensure that it is a numerical value
  # and that the name of the variable is indeed 'expr'
    if (!is.double(d$expr)) {
    stop("'expr' must be the name of an existing numerical value.")
  }
  

  # Parser technical variation, ensure that it is a factor with only 2 levels
  # and that its name is 'type'
  if (!is.factor(d$type)) {
    stop("'type' must be the name of an existing factor, 
         whose first level refers to the pool-and-split controls.")
  }
  if (length(levels(d$type)) != 2) {
    stop("The technical variation must be a factor with strictly 2 levels, 
         the first one referring to the pool-and-split controls, 
         and the second one to the 10-cell samples.")
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
  d.mm <- data.frame(stats::model.matrix(~ ., data = dplyr::select(d.cv, -id))[,-1])
  # Keep only the relevant variables in a tidy format
  # From now on, 'D' is the result of arranging and checking the input
  # provided by the user
  D <- dplyr::select(d, id, gene, expr, type) %>%
    dplyr::left_join(dplyr::mutate(d.mm, id = 1:dplyr::n()))
  # 
  N.total.genes <- length(unique(D$gene))
    
  ###########################################       
  ## Run the chosen steps of the pipeline: 
  ##   i) screening
  ##  ii) testing
  ## iii) correcting
  ##  iv) inferring
  ###########################################       
  
  ###########################################       
  ## Step i: screening
  ## 
  ## Estimate a GAMLSS;
  ## keep only those genes where samples have larger technical variation than controls;
  ## and where the p-value (which p-value?!!) is <= 0.05 after Sidak correction
  ## 
  ###########################################
      
if (screening) {                  
  #-----------------------------------------------------
  # First call to 'stochprofIC.gamlss()'
  # do it gene by gene in a loop, and put it all back together into
  # the object called 'out.gamlss1',
  # that will contain a structured dataframe with the minimal required results
  # (gene, sigma coefficient for technical variation, p-value for technical variation)
  #
  out.gamlss1 <- NULL
  list.genes <- levels(D$gene)
  for (g in 1:length(list.genes)) {
    gene.now <- list.genes[g]
    d.gene <- D %>%
      filter(gene == gene.now)
    out.gamlss1 <- bind_rows(out.gamlss1, 
                   bind_cols(gene = gene.now,
                   stochprofIC.gamlss(d.gene,
                               re = re, re.placement = re.placement,
                               log = log.1st,
                               add.1 = add1.1st,
                               trim = trim.1st, 
                               zeros = zeros.1st, 
                               zero.value = zero.value.1st)))
  }

  
  
  #-----------------------------------------------------
  # Second part: Run a Sidak Correction of sigma-type p-values and keep only
  # sigma.type > 0 & Pval.sidak <= 0.05
  
  # Manual test for extracting sidak p-values
  #sidak(out.gamlss1$Pval)
  #extract.sidak(out.gamlss1$Pval)
  
print("Is Pval known at this point?")  
  out.sidak <- out.gamlss1 %>%
    dplyr::mutate(Pval.sidak = extract.sidak(Pval))
  
  out.screening1 <- out.gamlss1 %>%
    dplyr::mutate(Pval.sidak = extract.sidak(Pval)) %>%
    dplyr::filter(Estimate >= 0) %>%
    dplyr::filter(Pval.sidak <= 0.05)
  
  # Proportion of genes kept after sidak correction and removing negative sigma estimates
  message(paste("Proportion of genes kept after sidak correction\nand removing negative sigma estimates:", (dim(out.screening1)[1] / N.total.genes)*100, "%"))
}

  
  ###########################################       
  ## Step ii: testing
  ## 
  ## Estimate a GAMLSS to identify possible batch effects
  ## 
  ###########################################

if (testing) {
  
  #-----------------------------------------------------
  # Third part: Run a GAMLSS with log of expr+1.
  # Return a single row for each gene with the Sigma coefficient and sigma p-value
  # for the variable that captures the technical/biological variation (type) as well as
  # the rest of mu and sigma coefficients to correct for technical variation
  # It requires 
  
  # Manual test for a single gene
  #stochprofIC.gamlss(filter(d, gene == levels(d$gene)[1]))
  
  # Run several genes
  d.post.sidak <- D %>%
    dplyr::filter(gene %in% out.screening1$gene) %>%
    dplyr::mutate(gene = factor(gene)) %>% 
    droplevels()
  out.gamlss2 <- NULL
  for(g in 1:length(levels(d.post.sidak$gene))){
    gene.now <- levels(d.post.sidak$gene)[g]
    out.gamlss2 <- dplyr::bind_rows(out.gamlss2, 
          dplyr::bind_cols(gene = levels(d.post.sidak$gene)[g],
            stochprofIC.gamlss(dplyr::filter(d.post.sidak, 
                       gene == levels(d.post.sidak$gene)[g]),
                       re = re, re.placement = re.placement,
                       log = log.2nd,
                       add.1 = add1.2nd,
                       trim = trim.2nd, 
                       zeros = zeros.2nd, 
                       zero.value = zero.value.2nd,
                       return = "simple")))
  }
  
  #-----------------------------------------------------
  # Fourth part: Run a FDR Correction of sigma-type p-values and keep only
  # sigma.type > 0 & Pval.fdr <= 0.05
  
  # Manual test for extracting bonferroni p-values (FDR)
  out.fdr <- out.gamlss2 %>%
    dplyr::mutate(Pval.dfr = stats::p.adjust(Pval, method="fdr")) 
  out.fdr
  
  out.screening2 <- out.gamlss2 %>%
    dplyr::mutate(Pval.fdr = stats::p.adjust(Pval, method="fdr")) %>%
    dplyr::filter(Estimate >= 0) %>%
    dplyr::filter(Pval.fdr <= 0.05)
  out.screening2
  
  # Proportion of genes kept after FDR correction and removing negative sigma estimates
  message(paste("Proportion of genes kept after FDR correction\nand removing negative sigma estimates:", (dim(out.screening2)[1] / N.total.genes)*100, "%"))
  
}

###########################################       
## Step iii: correcting
## 
## Use the estimated GAMLSS coefficients from the previous testing step
## for batch correction
## (this only makes sense if testing was done)
## 
###########################################

if (correcting) {
  
  
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
    out.full <- stochprofIC.gamlss(d.now,
                           return = "full", 
                           re = re, re.placement = re.placement,
                           log = log.2nd,
                           add.1 = add1.2nd,
                           trim = trim.2nd, 
                           zeros = zeros.2nd, 
                           zero.value = zero.value.2nd)

    out.corrected <- dplyr::bind_rows(out.corrected, 
                               tibble::tibble(gene = levels(d.post.fdr$gene)[g],
                                      out.full$exprs))   # Christiane: maybe change the name "exprs" in the gamlss function -- a bit confusing
    out.coefficients <- dplyr::bind_rows(out.coefficients, 
                                  tibble::tibble(gene = levels(d.post.fdr$gene)[g],
                                  out.full$coefficients))
  }
  
  out.corrected <- out.corrected %>%
    dplyr::mutate(gene = factor(gene))
  out.coefficients <- out.coefficients %>%
    dplyr::mutate(gene = factor(gene),
           Variable = factor(Variable))
  
# Columns have to be renamed again!  
  
  
  return(list(Coefficients = out.coefficients,
              `Corrected output` = out.corrected))

}
}

