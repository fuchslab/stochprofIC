# make all required functions known

source("R/stochprofIC.gamlss.R")
source("R/extract.sidak.R")
source("R/extract.bonferroni.R")
source("R/stochprofIC.pipeline.R")

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

stochprofIC.pipeline(d = breast.cancer,
                TPM = "TPM",
                technical.variation = "type",
                covariates = "biopsy",
                gene = "gene", 
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
                screening=T,testing=T,correcting=T,inferring=T)

### If we want to save the full output of the stochprofIC.pipeline
full.output <- stochprofIC.pipeline(d = breast.cancer,
                                    TPM = "TPM",
                                    technical.variation = "type",
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
               TPM = "TPM",
               technical.variation = "type",
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


  
