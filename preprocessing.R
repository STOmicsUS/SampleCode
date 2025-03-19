library(Matrix)
library(Seurat)
library(dplyr)
library(tidyverse)
library(purrr)
library(magrittr)
library(patchwork)
library(stringr)
library(SingleCellExperiment)  
library(SingleR)
library(pheatmap)
library(biomaRt)

set.seed(123)
options(future.globals.maxSize = 3 * 1024^3)

sample_ID <- c("sample1", "sample2")
sample_Name <- c("sample1", "sample2")
file_dir <- "./V1.3/"

binsize <- c("200","100","50","20","cell")

file_paths<-list()
file_names<-list()
for (i in seq_along(sample_ID)) {
  for (j in seq_along(binsize)) {
    file_path <- paste0(file_dir,sample_ID[i], "/Y01041NE_",sample_Name[i],"_bin", binsize[j], "_addimage/addimage/Y01041NE.addimg.rds")
    file_paths <- c(file_paths, file_path)
    file_name <-  paste0("Y01041NE_",sample_Name[i],"_bin", binsize[j])
    file_names <- c(file_names, file_name)
  }
}


preProcessing <- function(rds_data,nMAD){
  metadata<-rds_data@meta.data
  median_nCount <- median(metadata$nCount_Spatial, na.rm = TRUE)
  abs_dev_nCount <- abs(metadata$nCount_Spatial - median_nCount)
  mad_nCount <- median(abs_dev_nCount, na.rm = TRUE)
  
  median_nFeature <- median(metadata$nFeature_Spatial, na.rm = TRUE)
  abs_dev_nFeature <- abs(metadata$nFeature_Spatial - median_nFeature)
  mad_nFeature <- median(abs_dev_nFeature, na.rm = TRUE)
  
  lower_nCount <- median_nCount - nMAD * mad_nCount
  upper_nCount <- median_nCount + nMAD * mad_nCount
  
  lower_nFeature <- median_nFeature - nMAD * mad_nFeature
  upper_nFeature <- median_nFeature + nMAD * mad_nFeature
  
  rds_data.filter <- subset(rds_data, subset = nCount_Spatial  > lower_nCount &
                              nCount_Spatial < upper_nCount &
                              nFeature_Spatial  > lower_nFeature &
                              nFeature_Spatial < upper_nFeature &
                              percent.mito < 10)
  
  return (rds_data.filter)
  
}

for (i in seq_along(file_paths)) {
  file_path <- file_paths[[i]]
  file_name <- file_names[[i]]
  print(paste0("++++++++++Processing " ,file_path))
  print(paste0("++++++++++reading ", file_name))
  st_data <- readRDS(file_path)
  print(paste0("++++++++++pre-processing ", file_name))
  st_data.filter <- preProcessing(st_data,5)
  print(paste0("++++++++++SCT ", file_name))
  st_data.norm <- Seurat::SCTransform(st_data.filter,assay="Spatial", vars.to.regress = "percent.mito", return.only.var.genes = FALSE,method = "glmGamPoi")
  print(paste0("++++++++++Saving ", file_name))
  saveRDS(st_data.norm, file = paste0(file_name,"SCT.rds"))
}


