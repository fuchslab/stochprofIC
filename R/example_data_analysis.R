# make all required functions known

source("R/stochprofIC.gamlss.R")
source("R/extract.sidak.R")
source("R/extract.bonferroni.R")
source("R/stochprofIC.pipeline.R")

require(gamlss)
require(dplyr)

# load breast cancer data and apply pipeline

load("data/breast-cancer-data-ready-to-pipeline.RData")
View(breast.cancer)

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
                re.placement = "both")


# load glioma data and apply pipeline
# 
load("data/glioma-data-ready-to-pipeline.RData")
View(glioma)

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
  
