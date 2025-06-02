################################
### Prepare the environement ###
################################

# Libraries
library(loupeR)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(dplyr)

# Functions and palettes
source(file = "/scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/my_R_functions/Dolgos_2024_Custom_Functions.R")

# Create output path
out_path = create_exp_folder(
  github_dir = github_dir,
  samples_ID = "Dolgos_2024_Early_Passage_Henry_Tabula",
  exp = "5_final_export"
)

#####################
### Load the data ###
#####################

# Load sce after normalization (no rescaling)
file_path = get_exp_file_path(
  github_dir = github_dir,
  samples_ID = "Dolgos_2024_Early_Passage_Henry_Tabula",
  prev_exp = "3_dimension_reduction", 
  pattern = "sce_comb_dim_reduced.rds" # Sce after batch correction for visualization
)

sce_comb = readRDS(file = file_path)

##############################################
### Select and clean the exported metadata ###
##############################################

# Ensure the UMAP coordinates are stored in colData of the SCE object
coords <- as_tibble(reducedDim(sce_comb, "UMAP_on_MNN.1-20"))[, 1:2]
colData(sce_comb)$UMAP1 <- coords[[1]]  # Store the first column as UMAP1
colData(sce_comb)$UMAP2 <- coords[[2]]  # Store the second column as UMAP2


# Extract the desired colData and convert to a tibble
colData_table <- as.data.frame(colData(sce_comb))

# Export per cell sce metadata
colnames(colData_table) = colnames(colData(sce_comb))
rownames(colData_table) = colnames(sce_comb)

# Export the colData table
write.csv(colData_table, 
          file = paste0(out_path, "Dolgos_2024_Early_Passage_Henry_Tabula_metadata.csv"),
          quote = FALSE, row.names = TRUE, col.names = TRUE
)

# Export cleaned sce_comb object
saveRDS(
  object = sce_comb,
  file = paste0(out_path,"Dolgos_2024_Early_Passage_Henry_Tabula_sce_final.rds"))

###################################################
### Export matrix, features, and barcode files ###
###################################################

# Extract the count matrix
counts_matrix <- counts(sce_comb)  
counts_matrix <- Matrix(counts_matrix, sparse = TRUE)

# Export the count matrix
writeMM(counts_matrix, 
        file = paste0(out_path, "Dolgos_2024_Early_Passage_Henry_Tabula_matrix.mtx")
)

# Export features.tsv (gene names)
write.table(rowData(sce_comb)$SYMBOL, 
            file = paste0(out_path, "Dolgos_2024_Early_Passage_Henry_Tabula_features.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Export barcodes.tsv (cell barcodes)
write.table(colnames(sce_comb), 
            file = paste0(out_path, "Dolgos_2024_Early_Passage_Henry_Tabula_barcodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)


###########################
### Create a loupe file ###
###########################

# Convert sce to seurat object
seurat_comb = as.Seurat(
  x = sce_comb, 
  counts = "counts", 
  data = "logcounts"
)

# Export Seurat 
saveRDS(
  object = seurat_comb,
  file  = paste0(out_path,"Dolgos_2024_Early_Passage_Henry_Tabula_seurat_final.rds")
)

# Fix colnames suffix repalcing "-" by "_" to work with LoupeR
colnames(seurat_comb) = gsub(x = colnames(seurat_comb), pattern = "-", replacement = "_")

# Extract count matrix, clusters and projections for loupe file
assay = seurat_comb[["originalexp"]]

# Select only useful clusters present in the seurat object
clusters = select_clusters(seurat_comb)

# Select only useful projections present in the seurat object
projections = loupeR::select_projections(seurat_comb)
projections = projections[which(names(projections) == c("UMAP_on_MNN.1.20"))]


# Create loupe file
create_loupe(
  assay@counts,
  clusters = clusters,
  projections = projections,
  output_dir = out_path,
  output_name = "Dolgos_2024_Early_Passage_Henry_Tabula_loupe",
  force = T)

