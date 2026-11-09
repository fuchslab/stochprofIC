## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
library(stochprofIC)

## -----------------------------------------------------------------------------
data(breast.cancer)
breast.cancer

## ----echo = TRUE, results = "hide", warning = FALSE, message = FALSE, error = FALSE----
screening.output <- stochprofIC.screening(
                d = breast.cancer, # Working dataset
                expr = "TPM", # Name of the column containing the expression variable
                gene = "gene",  # Name of the column containing the gene variable
                type = "type", # Name of the column containing the type variable
                covariates = "biopsy", # Name(s) of the column(s) containing the covariates to be included
                #re = NA,   # not implemented yet, refers to Random Effects
                log = FALSE, # Whether to log the expression in the screening gamlss call
                add1 = FALSE,  # Whether to add +1 to the expression in the screening gamlss call
                trim = 3, # Number of values to be trimmed in the screening gamlss call
                zeros = TRUE, # Whether to change the value of zero in expr in the screening gamlss call
                zero.value = NA#, # Value to substitute the zero in expr in the screening gamlss call
                #re.placement = "both") # Not yet implemented, refers to Random Effects
)

## -----------------------------------------------------------------------------
screening.output

## ----echo = TRUE, results = "hide", warning = FALSE, message = FALSE, error = FALSE----
testing.output <- stochprofIC.testing(d = breast.cancer, # Working dataset
                                       expr = "TPM", # Name of the column containing the expression variable
                                       gene = "gene", # Name of the column containing the gene variable
                                       type = "type", # Name of the column containing the type variable
                                       covariates = "biopsy", # Name(s) of the column(s) containing the covariates to be included
                                       re = NA,   # not implemented yet, refers to Random Effects
                                       log = TRUE, # Wether to log the expression in the testing gamlss call
                                       add1 = TRUE, # Wether to add +1 to the expression in the testing gamlss call
                                       trim = NA, # Number of values to be trimmed in the testing gamlss call
                                       zeros = FALSE, # Whether to change the value of zero in the testing gamlss call
                                       zero.value = 0, # Whether to substitute the zero in the expr in the testing screening gamlss call
                                       re.placement = "both") # Not yet implemented, refers to Random Effects

## -----------------------------------------------------------------------------
testing.output

## ----echo = TRUE, results = "hide", warning = FALSE, message = FALSE, error = FALSE----
correcting.inferring.output <- stochprofIC.pipeline(d = breast.cancer,
                                       expr = "TPM",
                                       gene = "gene", 
                                       type = "type",
                                       covariates = "biopsy",
                                       re = NA,   # not implemented yet
                                       log.1st = FALSE, 
                                       log.2nd = TRUE,
                                       add1.1st = FALSE, 
                                       add1.2nd = TRUE,
                                       trim.1st = 3, 
                                       trim.2nd = NA,
                                       zeros.1st = TRUE, zeros.2nd = FALSE,
                                       zero.value.1st = NA, zero.value.2nd = 0,
                                       re.placement = "both")

## -----------------------------------------------------------------------------
correcting.inferring.output

## -----------------------------------------------------------------------------
correcting.inferring.output$`Corrected output`
correcting.inferring.output$`Coefficients`

## ----echo = TRUE, results = "hide", warning = FALSE, message = FALSE, error = FALSE----
full.output <- stochprofIC.pipeline(d = breast.cancer,
                                    expr = "TPM",
                                    type = "type",
                                    covariates = "biopsy",
                                    gene = "gene")

## -----------------------------------------------------------------------------
full.output$Coefficients
full.output$Corrected
full.output$`Genes with error`

