###################################################
## stochprofIC Examples on Data Preparation
###################################################
##
## This file contains two examples on how to transform the original data in csv
## into a tidy long format data suitable for the functions in the package
#
##required packages
#require(magrittr)
#require(dplyr)
#require(tidyverse)
#
##----------------------------------- EXAMPLE BREAST CANCER
## Load original data for breast cancer
## Rename the variable "symbol" that contains the genes into "gene"
## Select only the "gene" variable and samples with RSEM data
#bc <- read.csv("data/GSE147356_2020_02_06_breastcancer_allcompiled.csv") %>%
#  rename(gene = symbol) %>%
#  mutate(gene = as.factor(gene)) %>%
#  dplyr::select(-c(ensgene, entrez, chr, start, end, strand, biotype, description))
#
## Create a tidy dataset in the long format using pivot_longer
## In the same pipe, do the following steps:
## - Create specific variables from sample names (e.g. biopsy and type of sample)
## - Remove observations containing mitochondrial (MT) genes and ERCCs
## - Normalize RSEM to TPM (for each sample, thus grouping by sample)
#breast.cancer <- bc %>%
#  pivot_longer(-gene, names_to = "sample", values_to = "RSEM") %>%
#  mutate(sample = as.factor(sample)) %>%
#  mutate(temp = sample) %>%
#  separate(col=temp, into=c("biopsy", "type"),
#           sep="\\.") %>%
#  mutate(type = str_replace(type, "[:digit:]+", "")) %>%
#  mutate(biopsy = as.factor(biopsy),
#         type = as.factor(type)) %>%
#  filter(!grepl('^mt-', gene)) %>%
#  filter(!grepl('^ERCC-', gene)) %>%
#  group_by(sample) %>%
#  mutate(TPM = (RSEM / sum(RSEM) * 1000000)) %>%
#  ungroup()
#
## Save tidy dataset in the long format for subsequent steps
#save(breast.cancer, file="data/breast-cancer-data-ready-to-pipeline.RData")
#
##----------------------------------- EXAMPLE GLIOMA
## Load original data for glioma
## Rename the variable "symbol" that contains the genes into "gene"
## Select only the "gene" variable and samples with RSEM data
#gl <- read.csv("data/GSE147360_2020_03_20_opc_allcompiled.csv") %>%
#  rename(gene = symbol) %>%
#  mutate(gene = as.factor(gene)) %>%
#  dplyr::select(-c(ensgene, entrez, chr, start, end, strand, biotype, description))
#
## Load metadata for glioma with info on Sex, Time and Type of sample
#meta <- read.csv("data/GSE147360-GPL19057_series_matrix.txt", skip=29, sep="\t") %>%
#  filter(X.Sample_title == "!Sample_characteristics_ch1") %>%
#  tibble() %>%
#  dplyr::select(-X.Sample_title) %>%
#  pivot_longer(everything(), names_to = "sample", values_to="metadata") %>%
#  filter(!str_detect(metadata, "^genotype")) %>%
#  filter(!str_detect(metadata, "^section")) %>%
#  mutate(sex = case_when(
#    metadata == "Sex: female" ~ "female", 
#    metadata == "Sex: male" ~ "male",
#    TRUE ~ NA_character_)) %>%
#  mutate(time = case_when(
#    metadata == "time: 12 dpi (days following tamoxifen treatment start)" ~ "dpi12",
#    metadata == "time: 90 dpi (days following tamoxifen treatment start)" ~ "dpi90",
#    TRUE ~ NA_character_)) %>%
#  mutate(type = case_when(
#    metadata == "sample type: 10-cell sample" ~ "sample",
#    metadata == "sample type: pool-and-split control" ~ "control",
#    TRUE ~ NA_character_)) %>%
#  dplyr::select(-metadata) %>%
#  pivot_longer(-sample, names_to="variable", values_to = "values") %>%
#  filter(!is.na(values)) %>%
#  pivot_wider(names_from = "variable", values_from = "values") %>%
#  mutate(id = str_sub(sample, 1, 8))
#
## Create a tidy dataset in the long format using pivot_longer
## In the same pipe, do the following steps:
## - Create specific variables from sample names (e.g. sex, time and type of sample)
## - Remove observations containing mitochondrial (MT) genes and ERCCs
## - Normalize RSEM to TPM (for each sample, thus grouping by sample)
#glioma <- gl %>%
#  pivot_longer(-gene, names_to = "id", values_to = "RSEM") %>%
#  left_join(meta) %>%
#  filter(!is.na(sample)) %>%
#  dplyr::select(sample, gene, time, type, sex, RSEM) %>%
#  mutate(sample = as.factor(sample),
#         time = as.factor(time),
#         type = as.factor(type),
#         sex = as.factor(sex)) %>%
#  filter(!grepl('^mt-', gene)) %>%
#  filter(!grepl('^ERCC-', gene)) %>%
#  group_by(sample) %>%
#  mutate(TPM = (RSEM / sum(RSEM) * 1000000)) %>%
#  ungroup()
#
## Save tidy dataset in the long format for subsequent steps
#save(glioma, file="data/glioma-data-ready-to-pipeline.RData")
