# stochprofIC

stochprofIC is a tool to run a pipeline based on a generalized additive model to detect/parameterize transcriptional heterogeneity of 10-cell data.

Creator:

- Mercè Garí


Contributors:

- Mercè Garí
- Christiane Fuchs
- Kevin A. Janes


To install or update, run:
```
devtools::install_github("fuchslab/stochprofIC")
```
---

There are two vignettes with examples on how to clean and arrange datasets for the pipeline and on how to proceed with the pipeline that can be viewed using:
```
vignette('data_preparation')
vignette('data_analysis')
```

We demonstrate the usage of the pipeline on two subsets (50 genes) of the following data sets:

- ```glioma```

  which was published along this paper:
  
  Sutcliffe MD, Galvao RP, Wang L, Kim J, Rosenfeld LK, Singh S, Zong H, Janes KA. Premalignant Oligodendrocyte Precursor Cells Stall in a Heterogeneous State of Replication Stress Prior to Gliomagenesis. Cancer Res. 2021 Apr 1;81(7):1868-1882. doi: 10.1158/0008-5472.CAN-20-1037. Epub 2021 Feb 2. PMID: 33531372; PMCID: PMC8137536.
  
  paper available here:

  https://pubmed.ncbi.nlm.nih.gov/33531372/
  
  data available here:
  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147360
  
  (accessed as of Apr 06, 2021)

The dataset, ready for the pipeline, can be accessed by:
```
data(glioma)
```
  
- ```breast.cancer```

  which was published along this paper:
  
  Singh S, Sutcliffe MD, Repich K, Atkins KA, Harvey JA, Janes KA. Pan-Cancer Drivers Are Recurrent Transcriptional Regulatory Heterogeneities in Early-Stage Luminal Breast Cancer. Cancer Res. 2021 Apr 1;81(7):1840-1852. doi: 10.1158/0008-5472.CAN-20-1034. Epub 2021 Feb 2. PMID: 33531373; PMCID: PMC8137565.
  
  paper available here:
  
  https://pubmed.ncbi.nlm.nih.gov/33531373/
  
  data available here:
  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147356
  
  (accessed as of Jan 20, 2021)


The dataset, ready for the pipeline, can be accessed by:
```
data(breast.cancer)
```
