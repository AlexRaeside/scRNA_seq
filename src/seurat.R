# scRNA-Seq integration 

#
#
#
#



##------ Setting the environment ------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)

##------- PBMC Donor A ---------------------------------------------------------

##---- Normalize and Scale -----------------------------------------------------

# load in filtered data 
donorA_data <- Read10X(data.dir = "data/pbmc10x/filtered_matrices_mex/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
donorA <- CreateSeuratObject(
    counts = donorA_data, 
    project = "pbmc_donorA",
    min.cells = 3,
    min.features = 200)

# normalize 
donorA <- NormalizeData(donorA)

# scale 
donorA <- ScaleData(donorA)

# find variable features

donorA <- FindVariableFeatures(
    donorA, 
    selection.method = "vst", 
    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(donorA), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(donorA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
donorA_var_features <- plot1 + plot2

# linear dimension reduction 
donorA <- RunPCA(donorA, features = VariableFeatures(object = donorA))
donorA_loadings <- VizDimLoadings(donorA, dims = 1:2, reduction = "pca")

# heat map of PCA1
# in comparisons to the PBMC3x most genes are the same but some 
# real differences like RPL13 being important here but not MALAT

donorA_heatmap <- DimHeatmap(donorA, dims = 1, cells = 500, balanced = TRUE)

# Elbow plot to see how many PCAs are relevant
# the first 10
donorA_elbow_plot <- ElbowPlot(donorA)

##---- Cluster Cells  ----------------------------------------------------------


donorA <- FindNeighbors(donorA, dims = 1:10)
donorA <- FindClusters(donorA, resolution = 0.5)
donorA <- RunUMAP(donorA, dims = 1:10)
donroA_umpa <- DimPlot(donorA, reduction = "umap", label = TRUE) +
    ggplot2::ggtitle("DonorA")

##---- Finding marker genes ----------------------------------------------------

donorA_markers <- FindAllMarkers(
    donorA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
donorA_marker_genes_df <- donorA_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

write_csv(donorA_marker_genes_df, file = "tables/donorA_marker_genes.csv")

# looking at the marker genes I can tell:
#
# cluster 0 is a Naive CD4+ T based on CCR7 expression
#
# cluster 1 is probably a memory CD4 T based on IL7 but AQP3 is 
# harder to understand though it has been expressed in lymphomas, a B cell
# cancer. Very interesting! 
#
# cluster 4 is probably CD8 based on GZMK and DUSP2 being more highly expressed
# CD8 cells then CD4 t cells. GZMK found in NK cells 
# 
# cluster 2 is maybe NK cells or CD8 cells based on GZMH expression
# but CCL5 expression is common in both NK and CD8 cells 
#
# cluster 3 could also be NK cell based GNLY expression. Is more likely to be 
# NK cells based on both GNLY and NKG7 
#
# both CD8 T cells and NK cells originate from a 
# common lymphoid progenitor
#
# cluster 5 is likely a B cell based on CD79A expression 
#
# cluster 6 is likely CD14 or CD33 based on S100A9 expression. Highly 
# likely to be CD14+ Mono due to LYZ expression 
#
# cluster 7 is likely based on FCGR3A+ Mono FCGR3A expression 
#
# cluster 8 probably platelets or DK. Going to say cluster 8 
# is DK cells and unfortunately the platelet cells have been incorrectly
# labelled as as cluster 6 CD14+ Mono probably due to LYZ expression 



##--- EDA on the 8 clusters ----------------------------------------------------


donorA_featureplot <- FeaturePlot(donorA,
    features = c(
    "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))

# take a closer look at platelets and DC cells
# 4 genes expressed in platelets
# cluster 8 is DK cells and unfortunately the platelets have been 
# placed in cluster 7 

platelet_dc_genes <- c("PPBP", "FCER1A", "CST3")

donorA_cd_featureplot <- FeaturePlot(
    donorA,
    features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))


# NK cells can be grouped 

cluster_2_genes <- c("GZMH", "CCL5", "FGFBP2", "GNLY", "GZMB")
cluster_3_genes <- c("CLIC3", "GNLY", "NKG7", "PRF1", "GZMB")
typical_tcell_genes

FeaturePlot(
    donorA,
    features = unique(
        c(cluster_2_genes, cluster_3_genes, "CD8A", "NCAM1", "FCGR3A", "GZMB")))



x <- c("CD3G", "FCGR3A", "CD247", "TCRD")

FeaturePlot(
    donorA,
    features = x)



# the CD8A is expressed on 40% on NK cells 
# appears to be the difference betwwn the cells in cluster 2 and 3 here

