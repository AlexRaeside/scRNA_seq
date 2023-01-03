# scRNA-Seq integration 

# Initial analysis of 
# Frozen PBMCs (Donor A), Cell Ranger 1.1.0 dataset
#
#



##------ Setting the environment ------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(ggplot2)

##------- PBMC Donor A ---------------------------------------------------------

##---- Normalize and Scale -----------------------------------------------------

# load in filtered data 
pbmc <- Read10X(data.dir = "data/filtered_matrices_mex/hg19/")

gene_names <- data.frame(gene_name = pbmc@Dimnames[[1]])
write_csv(gene_names, file = "tables/pbmc_gene_names.csv")

# Initialize the Seurat object with the filtered non-normalized data
pbmc <- CreateSeuratObject(
    counts = pbmc, 
    project = "pbmc_donorA",
    min.cells = 3,
    min.features = 200)


# There is quite a lot of MT genes
# when comparing PBMC to BMPC 

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Each gene normalized by counts out of total counts in cell
# Pseudocount of 1 added to zero counts genes

pbmc <- NormalizeData(pbmc)

# Linear scaling so that every gene has a zero across all cells
pbmc <- ScaleData(pbmc)

# find variation
pbmc <- FindVariableFeatures(
    pbmc, 
    selection.method = "vst", 
    nfeatures = 2000)

# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(pbmc), 20)

top20
#  [1] "S100A9"   "S100A8"   "LYZ"      "IGLL5"    "CST3"     "FCER1A"   
# "HLA-DRA"  "PTGDS"    "CD74"     "HLA-DPA1"
# All seem fine for PBMC cells

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top20[1:10], repel = TRUE)

ggsave(
    filename = "figures/pbmc_donorA_variable_features.png",
    plot = plot2
)



# linear dimension reduction 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc_loadings <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")



# Elbow plot to see how many PCAs are relevant to the 
# deviation 
# the first 10
pbmc_elbow_plot <- ElbowPlot(pbmc)





##---- Cluster Cells  ----------------------------------------------------------

# K param nearest neighbors
# default to k = 20

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc_umpa <- DimPlot(
    pbmc, reduction = "umap", label = TRUE) +
    ggplot2::ggtitle("PBMC Donor A")


pbmc_umpa

ggsave(
    filename = "figures/pbmc_donorA_umpa.png",
    plot = pbmc_umpa
)


##---- Finding marker genes ----------------------------------------------------

pbmc_markers <- FindAllMarkers(
    pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_marker_genes_df <- pbmc_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

write_csv(
    pbmc_marker_genes_df, 
    file = "tables/pbmc_marker_genes.csv")



# looking at the marker genes I can tell:
#
# cluster 0 is a Naive CD4+ T based on CCR7 expression
#
# cluster 1 is probably a memory CD4 T based on IL7 but AQP3 is 
# harder to understand though it has been expressed in lymphomas, a B cell
# cancer. Very interesting! 
#
# Cluster 2 
# CCL5 expression is common in both NK and CD8 cells 
# 
# cluster 3 could be NK or CD8
#
# cluster 4 is probably CD8 based on GZMK and DUSP2 being more highly expressed
# CD8 cells then CD4 t cells. GZMK found in NK cells 
# 
# cluster 5 B Cells due to CD79A
#
# cluster 6 CD14+ Mono based on LYZ
#
# cluster 7 is likely FCGR3A+ Mono due to FCGR3A expression 
#
# cluster 8 probably platelets or DK. Going to say cluster 8 
# is DK cells and unfortunately the platelet cells have been incorrectly
# labelled as as cluster 6 CD14+ Mono probably due to LYZ expression 

new.cluster.ids <- c(
    "Naive CD4+ T", "Memory CD4 T", "CD8 T", "CD8 T / NK", "NK",
    "B", "CD14+ Mono", "FCGR3A+ Mono", "DC")




##--- Investigate Cluster 8 ----------------------------------------------------


# take a closer look at platelets and DC cells 
# FCER1A, CST3 are expressed in DK cells 
# and PPBP is expressed in platelets
# placed in cluster 7 

platelet_dc_genes <- c("PPBP", "FCER1A", "CST3")

donorA_cd_featureplot <- FeaturePlot(
    pbmc,
    features = platelet_dc_genes)

donorA_cd_featureplot

# unfortunately platelets appear to be placed into cluster 6 
# with CD14+ Mono

##--- Investigate Cluster 2 ----------------------------------------------------

# cluster 2 is the most problematic as could be ether CD8 T-cells
# or NK cells. Both of cells share a common progenitor so this is not
# surprising

# Lets look at some genes known to be expressed with the different cell types
#
# CD8A is a marker for CD8+ T
# GNLY (NKG5) and NKG7 are markers for NK cells 
# high expression of CCL5 is found in cluster 2 and known to be 
# more highly expressed in CD8 T cells then NK or CD4 T cells
# The gene FCGR3A for CD16 is commonly expressed in NK cells.
# The gene CD3G for CD3 is more associated with T cells then NK cells
# but subgroups of NK cells with CD3G expressed have been shown 


donor_nk_cd8_features <- FeaturePlot(
    pbmc,
    features = c("CD8A", "GNLY", "NKG7", "CCL5", "FCGR3A", "CD3G"))

donor_nk_cd8_features

# I am going to classify cluster 2 as a NK cell due to expression of 
# NKG7, GNLY and FCGR3A which are shared between clusters 2 and 3
# Due to the presence of CD8 in cluster 2 not cluster 3 I will refer to
# the cluster as CD8+_NK and cluster 3 as CD8-_NK


##--- Label Clusters -----------------------------------------------------------




names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc_labled_plot <- DimPlot(
    pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() +
    ggplot2::ggtitle("PBMC DonorA")

pbmc_labled_plot

ggsave("figures/labelled_pbmc_umap.png", pbmc_labled_plot)




# to help me see how the annotation lines up with the clustering 

cell_cluster_df <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("cell", "cluster", "cluster_cell_n")
colnames(cell_cluster_df) <- x


for (i in new.cluster.ids) {
    print(i)
    
    cells <- WhichCells(pbmc, idents = i)
    df <- data.frame(cell = cells)
    df$cluster <- rep(i, nrow(df))
    df$cluster_cell_n <- rep(nrow(df), nrow(df))
    
    
    cell_cluster_df <- rbind(cell_cluster_df, df)
    
    
}


write_csv(cell_cluster_df, file = "tables/pbmc_cell_clusters.csv")

saveRDS(pbmc, file = "rData/pbmc.rds")
