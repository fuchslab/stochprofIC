# make all required functions known

source("R/stochprofIC.pipeline.R")
source("R/stochprofIC.gamlss.R")
source("R/extract.sidak.R")
source("R/extract.bonferroni.R")

#install.packages("BiocManager")
#BiocManager::install("multtest")

require(gamlss)
require(dplyr)
require(multtest)

# load breast cancer data and apply pipeline

load("data/breast-cancer-data-ready-to-pipeline.RData")
View(breast.cancer)

# Select only 50 genes for testing
breast.cancer <- breast.cancer %>%
  filter(gene %in% levels(gene)[1:50]) %>%
  droplevels()

## for testing
#colnames(breast.cancer)[5] <- c("hamster")

screening.test <- stochprofIC.pipeline(d = breast.cancer,
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
                re.placement = "both",
                return="screening")
# This is correctly returning the first part of the pipeline (screening)
# It returns the genes that passes sidak, as well as the sigma type estimates from gamlss
screening.test


testing.test <- stochprofIC.pipeline(d = breast.cancer,
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
                                       re.placement = "both",
                                       return="testing")

# This is correctly returning the second part of the pipeline (testing)
# It returns the genes that passes FDR, as well as the sigma type estimates from gamlss
testing.test


correcting.inferring.test <- stochprofIC.pipeline(d = breast.cancer,
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
                                       re.placement = "both",
                                       return="correcting&inferring")
# This is correctly returning the last part of the pipeline (correcting&inferring)
# It returns a list with two objects:
# All the gamlss estimates from the final selected genes
# The corrected output data together with the non-corrected (TPM) for comparison.
correcting.inferring.test
# This second output (Corrected output) is the one that can be run in stochprofML
correcting.inferring.test$`Corrected output`




### If we want to save the full output of the stochprofIC.pipeline
### 
### The following call does not work yet because columns have to be re-renamed again.
full.output <- stochprofIC.pipeline(d = breast.cancer,
                                    expr = "TPM",
                                    type = "type",
                                    covariates = "biopsy",
                                    gene = "gene")

full.output$Coefficients
full.output$Corrected





# load glioma data and apply pipeline
# 
load("data/glioma-data-ready-to-pipeline.RData")
View(glioma)

# Select only 50 genes for testing
glioma <- glioma %>%
  filter(time %in% "dpi12") %>%
  filter(gene %in% levels(gene)[1:50]) %>%
  droplevels()

stochprofIC.pipeline(d = glioma,
               expr = "TPM",
               type = "type",
               covariates = "sex",
               gene = "gene", 
               re = NA,
               log.1st = FALSE, 
               log.2nd = TRUE,
               add1.1st = FALSE, 
               add1.2nd = TRUE,
               trim.1st = 3, 
               trim.2nd = NA,
               zeros.1st = TRUE, zeros.2nd = FALSE,
               zero.value.1st = NA, zero.value.2nd = 0,
               re.placement = "both")


  
