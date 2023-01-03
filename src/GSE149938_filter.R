# scRNA-Seq integration 

# Analysis of 
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

##----- Load the counts table --------------------------------------------------


wblood <- read_csv("tables/transposed_GSE149938_umi_matrix.csv")

colnames(wblood) <- as.character(wblood[1, ])
wblood <- wblood[-1,]


# add back the gene names 
gene_names <- read_csv("data/gse149938/gene_names.csv")
rownames(wblood) <- colnames(gene_names)

rownames(wblood)

# unfortunately there is no way to process all the cells 
# try with 4000 

keep_cols <- sample(colnames(wblood), 2000)

wblood <- wblood %>%
    dplyr::select(all_of(keep_cols))
write_csv(wblood, "tables/subset_counts.csv")

rownames(wblood) <- colnames(gene_names)
# add to Seurat object 
wblood <- CreateSeuratObject(
    counts = wblood, 
    min.cells = 3, 
    min.genes = 200, 
    project = "GSE149938_scRNAseq")

wblood[["percent.mt"]] <- PercentageFeatureSet(wblood, pattern = "^MT-")
wblood <- subset(wblood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

saveRDS(wblood, "rData/filtered_wblood.rds")
