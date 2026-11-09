## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(stochprofIC)
library(tidyverse)

## -----------------------------------------------------------------------------
gl <- data.table::fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147360&format=file&file=GSE147360%5F2020%5F03%5F20%5Fopc%5Fallcompiled%2Ecsv%2Egz") %>%
  rename(gene = symbol) %>%
  mutate(gene = as.factor(gene)) %>%
  dplyr::select(-c(ensgene, entrez, chr, start, end, strand, biotype, description))

## -----------------------------------------------------------------------------
meta.glioma <- read.csv("data/GSE147360-GPL19057_series_matrix.txt",
                        skip=29, sep="\t") %>%
  filter(X.Sample_title == "!Sample_characteristics_ch1") %>%
  tibble() %>%
  dplyr::select(-X.Sample_title) %>%
  pivot_longer(everything(), names_to = "sample", 
               values_to="metadata") %>%
  filter(!str_detect(metadata, "^genotype")) %>%
  filter(!str_detect(metadata, "^section")) %>%
  mutate(sex = case_when(
    metadata == "Sex: female" ~ "female", 
    metadata == "Sex: male" ~ "male",
    TRUE ~ NA_character_)) %>%
  mutate(time = case_when(
    metadata == "time: 12 dpi (days following tamoxifen treatment start)" ~ "dpi12",
    metadata == "time: 90 dpi (days following tamoxifen treatment start)" ~ "dpi90",
    TRUE ~ NA_character_)) %>%
  mutate(type = case_when(
    metadata == "sample type: 10-cell sample" ~ "sample",
    metadata == "sample type: pool-and-split control" ~ "control",
    TRUE ~ NA_character_)) %>%
  dplyr::select(-metadata) %>%
  pivot_longer(-sample, names_to="variable", values_to = "values") %>%
  filter(!is.na(values)) %>%
  pivot_wider(names_from = "variable", values_from = "values") %>%
  mutate(id = str_sub(sample, 1, 8)) %>%
  filter(time %in% "dpi90") %>%
  dplyr::select(-time)

## ----message=FALSE------------------------------------------------------------
glioma <- gl %>%
  pivot_longer(-gene, names_to = "id", values_to = "RSEM") %>%
  mutate(id = str_replace(id, "-", ".")) %>%
  left_join(meta.glioma) %>%
  filter(!is.na(sample)) %>%
  dplyr::select(sample, gene, type, sex, RSEM) %>%
  mutate(sample = as.factor(sample),
         type = as.factor(type),
         sex = as.factor(sex)) %>%
  filter(!grepl('^mt-', gene)) %>%
  filter(!grepl('^ERCC-', gene)) %>%
  group_by(sample) %>%
  mutate(TPM = (RSEM / sum(RSEM) * 1e6)) %>%
  ungroup() %>%
  dplyr::select(sample, gene, type, sex, TPM)

## -----------------------------------------------------------------------------
glioma

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  #save(glioma, file="data/glioma-data-ready-to-pipeline.RData")

## -----------------------------------------------------------------------------
bc <- data.table::fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147356&format=file&file=GSE147356%5F2020%5F02%5F06%5Fbreastcancer%5Fallcompiled%2Ecsv%2Egz") %>%
  rename(gene = symbol) %>%
  mutate(gene = as.factor(gene)) %>%
  dplyr::select(-c(ensgene, entrez, chr, start, end, 
                   strand, biotype, description))

## -----------------------------------------------------------------------------
breast.cancer <- bc %>%
  pivot_longer(-gene, names_to = "sample", values_to = "RSEM") %>%
  mutate(sample = as.factor(sample)) %>%
  mutate(temp = sample) %>%
  separate_wider_delim(cols = temp, delim = ".", 
                       names=c("biopsy", "type")) %>%
  mutate(type = str_replace(type, "[:digit:]+", "")) %>%
  mutate(biopsy = as.factor(biopsy),
         type = as.factor(type)) %>%
  filter(!grepl('^mt-', gene)) %>%
  filter(!grepl('^ERCC-', gene)) %>%
  group_by(sample) %>%
  mutate(TPM = (RSEM / sum(RSEM) * 1e6)) %>%
  ungroup() %>%
  dplyr::select(sample, gene, type, biopsy, TPM)

## -----------------------------------------------------------------------------
breast.cancer

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  #save(breast.cancer, file="data/breast-cancer-data-ready-to-pipeline.RData")

