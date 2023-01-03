# scRNA-Seq integration 

# Proper anlaysis  
# GEO GSE149938
#
# Loading the data as a seurat object and filtering 



##------ Setting the environment ------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(Matrix)
library(tibble)

##----- Load the seurat object --------------------------------------------------


wblood <- readRDS("rData/filtered_wblood.rds")


##---- Normalize ---------------------------------------------------------------

wblood <- NormalizeData(wblood)
wblood <- FindVariableFeatures(
    wblood, 
    selection.method = "vst", 
    nfeatures = 2000)

all.genes <- rownames(wblood)
wblood <- ScaleData(wblood, features = all.genes)

#--- Find variable features ----------------------------------------------------

wblood <- RunPCA(wblood, features = VariableFeatures(object = wblood))

DimPlot(wblood, reduction = "pca")
ElbowPlot(wblood) # 12


wblood <- FindNeighbors(wblood, dims = 1:10)
wblood <- FindClusters(wblood, resolution = 0.5)

wblood <- RunUMAP(wblood, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(wblood, reduction = "umap")

saveRDS(wblood, file = "rData/wblood.rds")
