#prepare.data <- function() {
  
  # read the glioma data
  glioma <- read.csv(file="data/GSE147360_2020_03_20_opc_allcompiled.csv")
  
  # read the breast cancer data
  breastcancer <- read.csv(file="data/GSE147356_2020_02_06_breastcancer_allcompiled.csv")
  
#}