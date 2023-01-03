# Integrating PBMC and BMPC 

### ------- Set the environment -------------------------------------------------


library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

### ------- Integrate Seurat  -------------------------------------------------



seurat_list <- list(
    "pbmc" = readRDS("rData/pbmc.rds"),
    "bmbp" = readRDS("rData/bmpc.rds")
)

seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seurat_list)


immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

# Anchors and Seurat list no longer used
rm(immune.anchors)
rm(seurat_list)
gc()


# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)



p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
plots <- p1 + p2
ggsave("figures/integrated_umap.png", plots)


# add new metadata variable that splits PBMC and BMPC cells


sample <- case_when(
    immune.combined$orig.ident == "pbmc_donorA" ~ "PBMC",
    immune.combined$orig.ident != "pbmc_donorA" ~ "BMPC"
)

immune.combined$sample <- sample


p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
plots <- p1 + p2
ggsave("figures/sample_integrated_umap.png", plots)
ggsave("figures/cluster_integrated_umap.png", p2)

##---- Annotating Clusters ----------------------------------------------------


# To annotate clusters I will use both the previous annotation of the 
# bmpc cells and the expression of key marker genes 




# The cell_cluster_df will list which cells belong to which 
# clusters


cell_cluster_df <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("cell", "cluster", "cluster_cell_n")
colnames(cell_cluster_df) <- x

# Would usually use mappy 


for (i in seq(0,13)) {
    
    
    cells <- WhichCells(immune.combined, idents = i)
    df <- data.frame(cell = cells)
    df$cluster <- rep(i, nrow(df))
    df$cluster_cell_n <- rep(nrow(df), nrow(df))
    
    
    cell_cluster_df <- rbind(cell_cluster_df, df)
    
    
}


cell_cluster_df$sample <- case_when(
    str_detect(cell_cluster_df$cell, "bar") == TRUE ~ "BMPC",
    str_detect(cell_cluster_df$cell, "bar") == FALSE ~ "PBMC"
)


write_csv(cell_cluster_df, file = "pbmc_bmc_cell_clusters.csv")


bmpc_cells <- cell_cluster_df %>%
    dplyr::filter(sample == "BMPC")


bmpc_cells$type <- sapply(str_split(bmpc_cells$cell, "_"),"[[",1)


cluster_types_df <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("cluster", "types")
colnames(cluster_types_df) <- x


for (i in seq(0,13)) {
    print(i)
    
    sub_tbl <- bmpc_cells %>%
        dplyr::filter(cluster == i)
    
    types <- table(sub_tbl$type)
    
    df <- data.frame(
        "cluster" = i,
        "types" = types
    )
    
    cluster_types_df <- rbind(cluster_types_df, df)
}



write_csv(cluster_types_df, file = "pbmc_bmpc_cell_types.csv")


# looking at how the annotation matches up with clusters
# the results are less then ideal but is quite similar to
# how the BMPC cells clustered when analysed by themselves 

# Cluster 0 has 111 CD8T cell and 98 CD4T cells
# Cluster 1 has 91 CD4T cells
# Cluster 2 is mix of B cells 272 naiB, 142 regB cells 
# Cluster 3 is mostly CD8T
# Cluster 4 mostly kineNK
# Cluster 5 mostly CD8T
# Cluster 6 contains 126 NKP cells and 106 proB cells
# Cluster 7 NK cells
# Cluster 8 B cells, mostly preB but no PBMC cells. Could be interesting!
# Cluster 9 contains B cells mostly memB
# Cluster 10 contains a few Monoctyes
# Cluster 11 contains NK cells
# Cluster 12 plasma cells
# Cluster 13 contains CLP cells

common_pbmc_tcell_markers <- c(
    "IL7R", "CCR7", "S100A4", "CD8A")

tcel__marker_plot <- FeaturePlot(
    immune.combined, features = common_pbmc_tcell_markers)

ggsave("figures/tcel_markers_immune_combined.png", tcel__marker_plot)


# The expression of CCR7 and IL7R in Cluster 0 links the cluster to Naive CD4+ T
# The expression of IL7R and S100A4 in Cluster 1 links the cluster to Memory CD4+
# The expression of CD8A is sporadically expressed across clusters 5 and 3
# Likely both clusters are CD8+ T


gzmk__marker_plot <- FeaturePlot(
    immune.combined, features = c("GZMK", "GZMB"))

# More GZMK expression in the CD8+ Tcells in Cluster 5 then cluster 3
# Cluster 5 will be referred to as GZMK+ CD8+ Tcells
# Cluster 3 will be referred to as GZMK- CD8+ Tcells

# Identifying B Cell clusters

bcell_marker_plot <- FeaturePlot(
    immune.combined, features = c("MS4A1"))

bcell_marker_plot

# Looking at MS4A1 expression clusters 2 and 9 are clearly B cells
# Looking annotated cells in cluster 2 and 9 
# Cluster 9 will be referred to as Memory B cells 
# Cluster 2 will be referred to as Naive/Regulatory B cells
# Cluster 8 will be referred to as preB 


# Identifying NK cell clusters 
nk_markers <- c("GNLY", "NKG7")
nk_marker_plot <- FeaturePlot(
    immune.combined, features = nk_markers)

# There is NK expression in cluster 3 but more in 7 and 11
# Cluster 7 contains mostly kineNK and will labelled Cytokine NK
# Cluster 11 contains mostly toxiNK and will be labelled Cytotoxic NK


# Identifying Monocyte Cell clusters

monocyte_markers <- c("CD14", "LYZ", "FCGR3A", "MS4A7")
monocyte_marker_plot <- FeaturePlot(
    immune.combined, features = monocyte_markers)

# based on the marker expression
# cluster 10 is FCGR3A+ Mono
# cluster 4 is CD14+ Mono

# Identifying dendritic monocytes

dendritic_monocytes<- c("FCER1A", "CST3")
dendritic_monocytes_plot <- FeaturePlot(
    immune.combined, features = dendritic_monocytes)

# small amount of expression of FCER1A in cluster 4
# suggests dendritic cells now part of the CD14+ Mono cluster


# Other Clusters

# Cluster 6
# only 6 out of the 304 cells in cluster 6 are PBMC based samples
# the rest are from BMPC suggesting that the cluster contains 
# a bone marrow related cell type 
# NKP belongs to intermediate progenitor cluster
# Cluster 6 to be labelled to as NKP

# Cluster 12 
# only 2 out of the 39 cells in cluster 12 are from PBMC bases samples
# plasma cells are found PBMC samples though
# could be there was just a low amount of the original PBCM sample 
# was plasma cells. Eccentric nucleus doesn't really count as 
# mononuclear 

# Cluster 13
# Cluster 13 is another BMPC cluster 
# containing most CLP common lymphoid progenitor 



# Rename Idents


immune.combined <- RenameIdents(
    immune.combined, 
    `0` = "Naive CD4+ T", 
    `1` = "Memory CD4+", 
    `2` = "Naive/Regulatory B cells",
    `3` = "GZMK- CD8+ Tcells", 
    `4` = "CD14+ Mono",
    `5` = "GZMK+ CD8+ Tcells", 
    `6` = "NKP", 
    `7` = "Cytokine NK", 
    `8` = "preB",
    `9` = "Memory B cells ",
    `10` = "FCGR3A+ Mono",
    `11` = "Cytotoxic NK", 
    `12` = "Plasma cells", 
    `13` = "Common Lymphoid Progenitor")

p3 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("figures/labelled_integrated_umap.png", p3)

p4 <- DimPlot(
    immune.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = "sample")
ggsave("figures/sample_labelled_integrated_umap.png", p4)


p3
p4



## ----- Common features and differential expressed genes --------------------


# For each cluster with substantial proportion of cells in both PBMC and BMPC
# create a table of the differential expressed genes 

DefaultAssay(immune.combined) <- "RNA"

compare_cluster <- c(
    "Naive CD4+ T", 
    "Memory CD4+", 
    "Naive/Regulatory B cells",
    "GZMK- CD8+ Tcells", 
    "CD14+ Mono",
    "GZMK+ CD8+ Tcells", 
    "Cytokine NK", 
    "Memory B cells ",
    "FCGR3A+ Mono",
    "Cytotoxic NK")

# Seurat workflow compares the average expression of the genes between
# the two types of cells (PBMC and BMPC) 
# will show broadly difference in gene expression in cell types 
# between samples but there no P values attached 
# I will save the tables showing gene expression changes for each 
# cell type and briefly discuss any interesting changes I find based  
# on the absoulte difference in expression

for(cluster in compare_cluster){
    
    print(cluster)
    
    
    subset.cells <- subset(immune.combined, idents = c(cluster))
    Idents(subset.cells) <- "sample"
    
    avg.cells <- as.data.frame(log1p(AverageExpression(subset.cells, verbose = FALSE)$RNA))
    avg.cells$gene <- rownames(avg.cells)
    avg.cells$difference <- abs(avg.cells$PBMC - avg.cells$BMPC)
    
    cluster <- str_replace_all(cluster, " ", "_")
    cluster <- str_replace_all(cluster, "/", "_")
    
    file_name <- paste0("tables/immune_combined_gene_expression_", cluster, ".csv")
    
    
    # There is a issue with different gene annotation between the two 
    # samples. PBMC is hg19 and BMPC is hg38
    # for most gene symbols this isn't a problem but for lncRNA 
    # and MT genes these change names quite often in updated genome
    # annotation 
    # the result is MALAT1 is differential expressed since it 
    # is present in PBMC genes and not in BMPC genes 
    # Ideally I would change the annotation in the BMPC dataset but 
    # as a quick fix I will only show gene which were expressed more 
    # then 0 in both samples 
    
    avg.cells <- avg.cells %>%
        dplyr::filter(PBMC != 0) %>%
        dplyr::filter(BMPC != 0)
    
    
    write_csv(avg.cells, file = file_name)
    
}








