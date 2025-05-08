
# Install
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

cran_packages <- c("dplyr", "ggplot2", "ggrepel", "igraph", 
  "openxlsx", "R.utils", "tidyverse", "msigdbr", "survminer",
  "survival", "WGCNA", "glmnet", "logistf", "randomForest",
  "pROC", "reticulate", "oncoPredict", "cowplot"
)
invisible(sapply(cran_packages, install_if_missing))



if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install ("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("impute")


# Load
library(dplyr)
library(ggplot2)
library(ggrepel)
library(igraph)
library(openxlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(R.utils)
library(tidyverse)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(msigdbr)
library(survminer)
library(survival)
library(WGCNA)
library(glmnet)
library(logistf)
library(randomForest)
library(pROC)
library(reticulate)
library(oncoPredict)
library(impute)
library(cowplot)