# scRNA-Seq integration 

# Initial analysis of 
# GEO GSE149938
#
#



##------ Setting the environment ------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(Matrix)
library(tibble)

##------- GSE149938 ------------------------------------------------------------

# As I was unable to find the filtered reads I will load the full table 
# and preform filtering 
# I may have to just use 5000 or the 7643 cells


raw_counts <- read_csv("data/gse149938/GSE149938_umi_matrix.csv")

head(colnames(raw_counts))
# [1] "OR4F5"      "FO538757.3" "FO538757.2" "OR4F29"     "OR4F16"     "SAMD11"  

head(rownames(raw_counts))
# [1] "1" "2" "3" "4" "5" "6"

head(raw_counts[,1])
# # A tibble: 6 × 1
#OR4F5             
#<chr>             
#1 BNK_spBM1_L1_bar25
#2 BNK_spBM1_L1_bar26
#3 BNK_spBM1_L1_bar27

# there are three problems in the data I have downloaded:
# 1. Gene names might be different to my PBMC dataset this is hg38 and the
# PBMC is hg17 
# 2. Data needs to be transposed with genes as rows 
# 3. The cell names have been loaded under the gene name OR4F5

# fix the col names 
cols <- c("cell", colnames(raw_counts))

# reload the raw counts skipping the headers 

raw_counts <- read_csv(
    "data/gse149938/GSE149938_umi_matrix.csv",
    col_names = FALSE,
    skip = 1)
colnames(raw_counts) <- cols


# write this as computer might crash 

write_csv(raw_counts, "tables/fixed_GSE149938_umi_matrix.csv")


# move col to rowname

raw_counts <- raw_counts %>%
    as.matrix() %>%
    t()

write_csv(as.data.frame(raw_counts), "tables/transposed_GSE149938_umi_matrix.csv")
