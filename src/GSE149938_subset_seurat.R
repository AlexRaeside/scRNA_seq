## Seurat analysis of subset_GSE149938_umi_matrix.csv


##------ Setting the environment ------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(ggplot2)



##------- BMPC -----------------------------------------------------------------


##------ Loading the data as a Seurat Object -----------------------------------

bmpc <- read_csv(
    "tables/subset_GSE149938_umi_matrix.csv",
    col_names = FALSE)

cells <- bmpc[,1]
bmpc <- bmpc[,-1]

# Will need genes as rows and cells as cols for Seurat

bmpc <- bmpc %>%
    t()

# add gene names as rows 
genes <- read_csv(col_names = FALSE, "tables/GSE149938_genes.csv")
genes <- unlist(genes[1,] %>% purrr::transpose())
rownames(bmpc) <- genes

# cell names as cols 
colnames(bmpc) <- cells$X1

# convert to seurat object 
bmpc<- CreateSeuratObject(
    counts = bmpc, 
    min.cells = 3, 
    min.genes = 200, 
    project = "BMPC")

## ---------- Filter   ------------------------------------------------


# during cell death mitochondria gene expression increases

bmpc[["percent.mt"]] <- PercentageFeatureSet(bmpc, pattern = "^MT-")
bmpc <- subset(
    bmpc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

saveRDS(bmpc, "rData/bmpc.rds")

bmpc <- readRDS("rData/bmpc.rds")

##------ - Normalize -----------------------------------------------------------

bmpc <- NormalizeData(bmpc)

##------- FindVariableFeature --------------------------------------------------

bmpc <- FindVariableFeatures(bmpc, selection.method = "vst", nfeatures = 2000)


## ------ Scale  ---------------------------------------------------------------


all.genes <- rownames(bmpc)
bmpc <- ScaleData(bmpc, features = all.genes)

##----- Get Clusters -----------------------------------------------------------


bmpc <- RunPCA(bmpc, features = VariableFeatures(object = bmpc))


# see some key genes like LYZ, S100A9, S100A8
# all fine


ElbowPlot(bmpc)
# PC 11 

bmpc <- FindNeighbors(bmpc, dims = 1:10)
bmpc <- FindClusters(bmpc, resolution = 0.5)


bmpc <- RunUMAP(bmpc, dims = 1:10)
bmpc_plot <- DimPlot(bmpc, reduction = "umap")

bmpc_markers <- FindAllMarkers(
    bmpc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
bmpc_marker_genes_df <- bmpc_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

write_csv(bmpc_marker_genes_df, file = "tables/bmpc_marker_genes.csv")
ggsave("figures/bmpc_umpa.png", bmpc_plot)


##---- Use preexisting annotation on clusters ----------------------------------

p1 <- DimPlot(
    bmpc, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(bmpc, reduction = "umap", label = TRUE, repel = TRUE)
plots <- p1 + p2

ggsave("figures/bmpc_umpa.png", plots)

WhichCells(bmpc, idents = 10)


##---- Nice table for cells and clusters --------------------------------------

# to help me see how the annotation lines up with the clustering 

cell_cluster_df <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("cell", "cluster", "cluster_cell_n")
colnames(cell_cluster_df) <- x


for (i in seq(0,10)) {
    print(i)
    
    cells <- WhichCells(bmpc, ident = i)
    df <- data.frame(cell = cells)
    df$cluster <- rep(i, nrow(df))
    df$cluster_cell_n <- rep(nrow(df), nrow(df))
    
    
    cell_cluster_df <- rbind(cell_cluster_df, df)
    
    
}


write_csv(cell_cluster_df, file = "bmpc_cell_clusters.csv")

saveRDS(bmpc, "rData/bmpc.rds")

