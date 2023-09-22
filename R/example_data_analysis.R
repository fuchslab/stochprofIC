load("data/glioma-data-ready-to-pipeline.RData")
View(glioma)

pipeline_global(d = glioma,
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
  