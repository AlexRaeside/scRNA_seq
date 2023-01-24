# scRNA_seq

### The tidy up bracnh 

Hello and welcome to the tidy up branch. I have a somewhat okay nice basic project here but it needs tidying up. 
How am i going to to do this? The project can be split into three parts, four if end up doing a Quarto report. 

First part is looking at the PBMC Donor A data. The code is straight forward Seurat pipeline. Loading 10x ranger
data and peforming manual annotation on single cell dataset. Im not intreasted in the QC as the 10x webiste has that fine.
What I need to do is tidy the code and output figures showing the expression of all the marker genes I used in manual annoatation. 



## Project outline

*Objective*

Using Frozen PBMCs (Donor A), Cell Ranger 1.1.0 dataset available at 10x Genomics, 
please perform its integration with the whole blood dataset from GEO GSE149938 using 
scanpy (preferably, but please use R if it is more convenient for you). 

Please organize all the commands and results in the form of Jupyter notebook. 
Use the available pre-calculated filtered expression matrices.
Briefly discuss the differences and similarities of the two datasets on gene and cell level.

*Conclusion*

The datasets were integrated correcty using Seurat so the labels established in the whole 
blood dataset could be used to annotate the cells the PBMC donor A. The nodualrity of the 
results were also improved with cell types in a sinlge cluster when PBMC was analysed independetly 



## Data

### Donor A PBMC Cell Ranger

The [Donor A PBMC Cell Ranger](https://www.10xgenomics.com/resources/datasets/frozen-pbm-cs-donor-a-1-standard-1-1-0) 
contains Peripheral Blood Mononuclear cells are the blood cells with a round nucleus
which includes T, B, NK and monocytes cells. Excludes erythrocytes (red blood
cells) which have no nuclei and the granulocyte which have multi-lobed nuclei.
The PBMC are extracted using a hydrophobic colloid and density gradient 
centrifuge. Plasma the top layer, PBMCs the middle layer and 
polynuclear cell the bottom layer. 
A typical composition of PBMC includes 70% T cells, 15% monocyte/macrophages,
10% B cells and 15% NK cells. Platelets are typically present in PBMC
samples, though they do play a role in the immune response they do not 
contain a nucleus. 

Looking at the 10x analysis page, the PBMC donor A sample looks very similar to 
other 10x PBMC sample in the Seurat tutorial.  The scRNA-Seq data contains 
3000 cells. 

Pre-filtered 10x data downloaded and placed in the data folder.

```

wget -0 data/frozen_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz \
https://cf.10xgenomics.com/samples/cell-exp/1.1.0/frozen_pbmc_donor_a/frozen_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz

```

### GSE149938 

GSE149938 is the Gene Expression Omnibus ID for a [large scale scRNA-Seq
analysis of 7,463 immunophenotypic blood cells across 21 healthy donors](https://academic.oup.com/nsr/article/8/3/nwaa180/5896476?login=false).
This is interesting as the almost twice as many cells as Donor A but 21 times 
as many donors. The paper is really interested in lncRNAs and the role 
they play in the cell differentiation, which is interesting as most
marker-based PBMC annotation I've looked at focuses on well known proteins 
like xxx and xxx for xxxx. This analysis of protein coding 
and lncRNA expression and transcription factor networking 
went into creating a cell type prediction model which 
can be found [online](http://scrna.sklehabc.com/). The prediction model 
takes a csv file as input with cells as columns and gene names as rows.
Would be cool to try to automatically annotate some of PBMC donor A cells.

The paper uses a R package abcCellmap which combines Seurat and 
[scmap](https://scmap.sanger.ac.uk/scmap/), 
made by Vladimir Kiselev at the Sanger.


Immunophenotypic cells includes mononuclear cells in Donor A PBMC 
as well as the Hematopoietic stem/progenitor cells (HSPCs). HSPCs
are the progenitor of both lymphoid and myeloid cells. This 
is because both cells from PBMC and harvested bone marrow were 
collected. Bone marrow contains 
CD34+ HSPCs, B cells, NK cells, T cells, monocytes, neutrophils
and erythrocytes. PBMC contains regulatory B, naive B, memory B, 
cytotoxic NK, cytokine NK and T cells. 

Since there should be no neutrophils, erythrocyte or HSPC cells in the PBMC 
I will only download the use these cells. The
names of the T cells, NK cell, B cells and monocytes can be downloaded with 
the following commands.

```
# B cell
wget -O data/GSM4793029_B.txt.gz \
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4793nnn/GSM4793029/suppl/GSM4793029%5FB%2Etxt%2Egz

# T Cell
wget -O data/GSM4793031_T.txt.gz \
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4793nnn/GSM4793031/suppl/GSM4793031%5FT%2Etxt%2Egz

# NK Cells
wget -O data/GSM4793030_NK.txt.gz \
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4793nnn/GSM4793030/suppl/GSM4793030%5FNK%2Etxt%2Egz


# Monoctyes 
wget -O data/GSM4793032_Mo.txt.gz \
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4793nnn/GSM4793032/suppl/GSM4793032%5FMo%2Etxt%2Egz

gunzip data/*.txt.gz

```

The count matrix contains 7643 cells as rows and 19812 genes as cols  

```

# download the counts 


# extract gene names and add the string sample to make it sure its read
# correctly 


# gene names are symbols the same as Donor A but Donor A is hg19 and 
# GSE149938 is hg32 so there may be differences 

# extract cell names 
awk -F , '{print $1}' data/GSE149938_umi_matrix.csv > tables/GSE149938_cells.csv


# combine lists of cells of intereast 

cat data/GSM47930* >> tables/cells_of_interest.txt

# a weird unicode symbol on the cells_of_interest.txt
# means the R script src/fix_cell_list.R was used to load and write out
# table. Intreastingly there are 9384 cell names in cell_list.txt and 
# there are only 7643 cells in the study. Some cells have two names as they 
# exist as parts of two clusters or annotation groups. Once filtered
# tables/subset_GSE149938_umi_matrix.csv contains 2965 cells. Slightly 
# less then half of the cells which makes sense looking at the figures.
# the src/fix_cell_list.R also creates a single line csv file 
# tables/cell_subset_GSE149938_umi_matrix.csv which the filtered gene
# counts are added to


# filter cells

grep -F -f tables/edited_cells_of_interest.txt data/GSE149938_umi_matrix.csv > \
tables/subset_GSE149938_umi_matrix.csv


shuf -n 500 tables/subset_GSE149938_umi_matrix.csv > tables/sampled_subsetGSE149938_umi_matrix.csv


cat tables/sampled_subsetGSE149938_umi_matrix.csv >> tables/cell_subset_GSE149938_umi_matrix.csv



```


## Differential expressed genes

Looking at the differential expressed genes and marking some interesting 
changes.

Unfortunately genome annotation issue have effected results of the 
differential expression analysis 


CD14+ Mono
PBMC	BMPC	gene	difference
0	4.75573081481199	DEFA3	4.75573081481199
4.68972869624599	0	MALAT1	4.68972869624599
0	4.45796807543227	LTF	4.45796807543227
4.60506807878426	0.239142906500432	FTH1	4.36592517228382
3.875024555039	0	MT-CO1	3.875024555039
3.64851714961893	0	MT-CO3	3.64851714961893

MALAT and mito genes expressed in PBMC

FGCR3A+ Mono same thing
PBMC	BMPC	gene	difference
5.23070746500656	0	MALAT1	5.23070746500656
4.13871751591105	0	MT-CO1	4.13871751591105
5.0360257029716	1.12012322932071	FTH1	3.91590247365089
3.70419854799964	0	MT-CO3	3.70419854799964
3.441687191202	0	MT-ND1	3.441687191202
3.80748201893328	0.397288232701891	RPL21	3.41019378623139
3.27968763970721	0	MT-CO2	3.27968763970721
3.16751203836744	0	MT-ND2	3.16751203836744
3.08820988516419	0	MT-CYB	3.08820988516419
2.92918295912805	0	MT-ND4	2.92918295912805
2.69346946429975	0	RP11-290F20.3	2.69346946429975
2.65647495762294	0	MT-ATP6	2.65647495762294

Cytokine NK
PBMC	BMPC	gene	difference
5.69792855262776	0	MALAT1	5.69792855262776
4.10262797747209	0	MT-CO1	4.10262797747209
3.75335866596993	0	MT-CO3	3.75335866596993
3.43445220981025	0	MT-CO2	3.43445220981025
3.34051283221371	0	MT-ND4	3.34051283221371
3.23573716777162	0	MT-CYB	3.23573716777162
3.14548925657424	0	MT-ND1	3.14548925657424
3.11488581116175	0	MT-ND2	3.11488581116175


Cytotoix NK
PBMC	BMPC	gene	difference
5.34917416196294	0	MALAT1	5.34917416196294
3.87570041656866	0	MT-CO1	3.87570041656866
3.77083407684561	0	MT-CO3	3.77083407684561
3.55916817555595	0	MT-ND1	3.55916817555595
3.23621892269994	0	MT-ND2	3.23621892269994
3.21924351344834	0	MT-CO2	3.21924351344834
3.08736159820905	0	MT-ND4	3.08736159820905


MALAT is a lncRNA that has been linked to vitamin deficiency and cancer 
Coud 



