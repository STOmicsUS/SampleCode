library(Matrix)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(sctransform)
library(reticulate)
library(SingleCellExperiment)
library(SingleR)
library(pheatmap)
library(biomaRt)
library(tidyverse)
library(harmony)
library(patchwork)
library(cowplot)
library(RColorBrewer)

set.seed(123)
options(future.globals.maxSize = 3 * 1024^3)


data_path ="./"
out_path = "./"
# 
sample_Name <- c("230_HG_1", "1084_HG_2")

st_data = c(paste0(data_path,"Y01041NE_",sample_Name[1],"_bin20SCT.rds"),
            paste0(data_path,"Y01041NE_",sample_Name[2],"_bin20SCT.rds"))


st_data_obj=c(readRDS(st_data[1]),readRDS(st_data[2]))

titles=c("230_HG_1","1084_HG_2")

find_clusters=function(norm_data,dimension_pca,k){
  norm_data <- RunPCA(norm_data, verbose = TRUE)
  norm_data <- FindNeighbors(norm_data, reduction = "pca", dims = 1:dimension_pca,k.param = k)
  norm_data <- RunUMAP(norm_data, reduction = "pca", dims = 1:dimension_pca)
  norm_data <- FindClusters(norm_data,algorithm=4,method = "igraph",cluster.name="leiden",resolution=0.8, verbosse = TRUE)
  return(norm_data)
  }


#Process first sample

st_data.norm=find_clusters(st_data_obj[[1]],20,20)
title=titles[2]
bin="Bin20"
ElbowPlot(st_data.norm)

de_markers <- FindAllMarkers(st_data.norm,only.pos=TRUE,logfc.threshold=0.5,min.pct = 0.1)
sorted_markers <- de_markers[order(de_markers$cluster, -de_markers$avg_log2FC), ]
top_genes <- sorted_markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) 
write.csv(top_genes,file=paste0(out_path,title,".",bin,".Top.markerGenes.csv"))
top_genes$gene <- factor(top_genes$gene, levels = unique(top_genes$gene))

#######Seurat Mapping

sc_data=readRDS("./sc.combined.harmony_fil_clean.rds")
scRNA_sce <- as.SingleCellExperiment(sc_data)
st_sce <- as.SingleCellExperiment(st_data.norm)
new_sc_data=Seurat::SCTransform(sc_data, verbose = TRUE,return.only.var.genes = FALSE,method = "glmGamPoi")

markers_genes<-read.csv(file =paste0("./marker_genes_Broad_labels_ext_LFC2.csv"), header=TRUE)
markers_genes=unique(markers_genes$gene)
print(length(markers_genes))

anchors = FindTransferAnchors(reference=new_sc_data, query = st_data.norm,normalization.method = "SCT", features = markers_genes)
predictions.assay <- TransferData(anchorset = anchors, refdata = new_sc_data$Broad_labels_ext, prediction.assay=TRUE, weight.reduction=st_data.norm[["pca"]], dims=1:20)
non.mapping <- c()
for(i in 1:dim(predictions.assay)[1]){ if(sum(predictions.assay@data[i,])==0) non.mapping <- c(non.mapping, rownames(predictions.assay)[i])}
predictions.assay@misc$non.mapping <- non.mapping
predictions.assay@misc$mapping <- setdiff(levels(as.factor(new_sc_data$Broad_labels_ext)), non.mapping)
st_data.norm[["predictions"]] <- predictions.assay
