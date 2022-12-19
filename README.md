# scRNA_seq

I have been given the following task:

Using Frozen PBMCs (Donor A), Cell Ranger 1.1.0 dataset available at 10x Genomics, please perform its integration with the whole blood dataset from GEO GSE149938 using scanpy (preferably, but please use R if it is more convenient for you). Please organize all the commands and results in the form of Jupyter notebook. Use the available pre-calculated filtered expression matrices. Briefly discuss the differences and similarities of the two datasets on gene and cell level.

What does this mean? And how am I going to do this? First lets break down the main buzz word.

## Intergation

A description of scRNA-Seq intergtation from the  

Here is a little on intergration taken from the lab that makes the R scRNA-seq analysis program [Seurat](https://satijalab.org/seurat/articles/integration_introduction.html). The joint analysis of two or more single-cell datasets poses unique challenges. In particular, identifying cell populations that are present across multiple datasets can be problematic under standard workflows. Seurat v4 includes a set of methods to match (or ‘align’) shared cell populations across datasets. These methods first identify cross-dataset pairs of cells that are in a matched biological state (‘anchors’), can be used both to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions. Though I am not sure what experimental conditions we would be looking at in this excercise? The PBMC dataset is from single Donor were the GSE149938 is from 21 healthy donors so there might be something there. It depends really what metadata we can get on the 21 healthy donors and donor A.

According to a [HBC training website](https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html) its important too look at the clustering without integration before deciding whether we need to perform any alignment. Do not just always perform integration because you think there might be differences - explore the data. If we had performed the normalization on both conditions together in a Seurat object and visualized the similarity between cells, we would have seen condition-specific clustering.

So in other things the first steps in this anlaysis will be to perform normlization and clustering on the two datasets. 

## A little about the data 

The PBMC data comes from [Zheng et al, “Massively parallel digital transcriptional profiling of single cells”](https://www.nature.com/articles/ncomms14049)
and can be downloaded from [10x](https://www.10xgenomics.com/resources/datasets/frozen-pbm-cs-donor-a-1-standard-1-1-0). PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells. In humans, the frequencies of these populations vary across individuals, but typically, lymphocytes are in the range of 70–90 %, monocytes from 10 to 20 %, while dendritic cells are rare, accounting for only 1–2 %. Analysis of PBMC is often performed with single-cell B-cell receptor sequencing (scBCR-seq) and single-cell T-cell receptor sequencing (scTCR-seq) based on the scRNA-seq libraries.
Analysis of the scBCR-seq and scTCR-seq could be done with [immcantation](https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html).

The GSE149938 data is a tad more complex. An experimental summary can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149938) and the data consists of the transcriptomes of 7,643 human blood cells covering 32 immunophenotypic cell types across 21 healthy donors. The orignal paper can be found [here](https://academic.oup.com/nsr/article/8/3/nwaa180/5896476?login=false). As always with GEO getting the file you want is a run around, at least for me.  




## What is scanpy
Scanpy is a scalable toolkit for analyzing single-cell gene expression data built jointly with anndata. It includes preprocessing, visualization, clustering, trajectory inference and differential expression testing. The Python-based implementation efficiently deals with datasets of more than one million cells.

It is like a python version of Suerat with AnnData being like python version of SumerizedExperiment. A way to keep count data across samples in the same place as experimental information about the samples. 

The tutorials break the scanpy pipeline into 5 parts:

### Clustering 

Going from counts to cell types using clustering and biomarker gene expression. This is the main part of the pipeline. 

### Visualization

Mostly making the types of figures seen in the Seurat clustering pipeline plus a dendrogram. No really intreast to learn matplotlib when I've spent so much time on ggplot2 in R so will mostly just go through this quick.

### Trajectory inference

According to the wiki 

Trajectory inference or pseudotemporal ordering is a computational technique used in single-cell transcriptomics to determine the pattern of a dynamic process experienced by cells and then arrange cells based on their progression through the process. Single-cell protocols have much higher levels of noise than bulk RNA-seq, so a common step in a single-cell transcriptomics workflow is the clustering of cells into subgroups. Clustering can contend with this inherent variation by combining the signal from many cells, while allowing for the identification of cell types. However, some differences in gene expression between cells are the result of dynamic processes such as the cell cycle, cell differentiation, or response to an external stimuli. Trajectory inference seeks to characterize such differences by placing cells along a continuous path that represents the evolution of the process rather than dividing cells into discrete clusters.

Not needed for integration so will come back to another time. 

### Integrating datasets

### Spatial data
Don't need to do this for the assignment so giving a miss for now but might tackle it with another datastet another time 


## Usefull Links

I have found a juptyer notebook by the Sanger center on the scanpy [pipeline](https://github.com/cellgeni/notebooks/blob/master/notebooks/new-10kPBMC-Scanpy.ipynb). It goes through filtering which I should be able to skip but the normlization, clustering and biomarker selection shoul be good.

The scanpy [tutorials](https://scanpy-tutorials.readthedocs.io/en/latest/index.html) are also pretty good.

