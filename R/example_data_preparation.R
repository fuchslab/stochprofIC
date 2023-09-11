##################################################
# stochprofIC Examples on Data Preparation
##################################################
#
# This file contains two examples on how to transform the original data in csv
# into a tidy long format data suitable for the functions in the package

#----------------------------------- EXAMPLE BREAST CANCER
# Load original data for breast cancer
# Rename the variable "symbol" that contains the genes into "gene"
# Select only the "gene" variable and samples with RSEM data
bc <- read.csv("GSE147356_2020_02_06_breastcancer_allcompiled.csv") %>%
  rename(gene = symbol) %>%
  mutate(gene = as.factor(gene)) %>%
  dplyr::select(-c(ensgene, entrez, chr, start, end, strand, biotype, description))

# Create a tidy dataset in the long format using pivot_longer
# In the same pipe, do the following steps:
# - Create specific variables from sample names (e.g. biopsy and type of sample)
# - Remove observations containing mitochondrial (MT) genes and ERCCs
# - Normalize RSEM to TPM (for each sample, thus grouping by sample)
breast.cancer <- bc %>%
  pivot_longer(-gene, names_to = "sample", values_to = "RSEM") %>%
  mutate(sample = as.factor(sample)) %>%
  mutate(temp = sample) %>%
  separate(col=temp, into=c("biopsy", "type"),
           sep="\\.") %>%
  mutate(type = str_replace(type, "[:digit:]+", "")) %>%
  mutate(biopsy = as.factor(biopsy),
         type = as.factor(type)) %>%
  filter(!grepl('^mt-', gene)) %>%
  filter(!grepl('^ERCC-', gene)) %>%
  group_by(sample) %>%
  mutate(TPM = (RSEM / sum(RSEM) * 1000000)) %>%
  ungroup()

# Save tidy dataset in the long format for subsequent steps
save(breast.cancer, file="breast-cancer-data-ready-to-pipeline.RData")

#----------------------------------- EXAMPLE GLIOMA
