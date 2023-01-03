# fix the cell list which has a weird unicode character 


library(readr)

data <- read_csv(col_names = FALSE, file = "tables/cells_of_interest.txt")

write.table(
    data, 
    file = "tables/edited_cells_of_interest.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)


genes <- read_csv(col_names = FALSE, "tables/GSE149938_genes.csv")
genes <- unlist(genes[1,] %>% purrr::transpose())
cell_genes <- c("cell_name", genes)

cell_genes <- paste(cell_genes, sep = ",", collapse = ",")

write_file(cell_genes, file = "tables/cell_names_genes.csv")

write_file(cell_genes, file = "tables/cell_subset_GSE149938_umi_matrix.csv")


