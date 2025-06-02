################################
### Prepare the environment ###
################################

# Set the Conda library path as the primary path
.libPaths(c("/scicore/home/wykopa75/parmen0000/miniconda3/envs/infercnv_ubuntu/lib/R/library", .libPaths()))

# Libraries
library(infercnv)
library(dplyr)
library(Seurat)

# Set the number of threads based on SLURM_CPUS_PER_TASK environment variable
num_threads <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

# Functions
source(file = "/scicore/home/wykopa75/GROUP/rparmentier/sc_RNAseq/my_R_functions/Dolgos_2024_Custom_Functions.R")

# Create output path
out_path = create_exp_folder(
  github_dir = github_dir,
  samples_ID = "Dolgos_PDOs-Dolgos_Tissue_Matching",
  exp = "infer_cnv"
)

#####################
### Load the data ###
#####################

# Load sce_comb after dimension reduction
file_path = get_exp_file_path(
  github_dir = github_dir,
  samples_ID = "Dolgos_PDOs-Dolgos_Tissue_Matching",
  prev_exp = "2_dimension_reduction", 
  pattern = "sce_comb_dim_reduced.rds")

print(file_path)
sce_comb = readRDS(file = file_path) 

# Convert sce to seurat object
seurat_comb = as.Seurat(
  x = sce_comb, 
  counts = "counts", 
  data = "logcounts"
)

#####################
### Run Infer CNV ###
#####################

seurat_comb$Unique_ID = paste0(seurat_comb$Sample_Description,"_", seurat_comb$Technology)

cat(unique(seurat_comb$Unique_ID))

infercnv_object <- CreateInfercnvObject(
  raw_counts_matrix = as.matrix(seurat_comb[["originalexp"]]@counts), # One should provide raw counts as infer_CNV performs it's own scaling operation
  annotations_file = as.matrix(seurat_comb$Unique_ID),
  gene_order_file = "~/cell_ranger_transcriptome_ref/infer_cnv_GRCh38_gene_ordering.txt",
  ref_group_names = c("P24-06_Org_Mat_10X_FRP"),
  chr_exclude = c("chrM"))


# Save reference object
saveRDS(
  object = infercnv_object,
  file = paste0(out_path, time_stamp(),"infercnv_object_reference_normal.rds")
)


infercnv_object <- infercnv::run(
  infercnv_object,
  num_threads = num_threads, # Number of threads (=cores) used to parellelize tasks
  cutoff = 0.1, # Better suits 10X data
  analysis_mode = "subclusters", # Doesn't perform subclustering on each sample
  num_ref_groups = 1, # We provided 6 type of references 
  cluster_references = F,
  min_cells_per_gene = 20,
  cluster_by_groups = T, # if TRUE Each group (= UniqueID) will be clustered separately, if FALS, follows k_obs_groups
  scale_data = T, # Better in the case of reference cell retrieved in an other context/technology
  HMM = T,
  denoise = T,
  plot_steps = T,
  useRaster = F, # If set T causes memory allocation problem with big datasets
  up_to_step = 100,
  out_dir= out_path)

# Save infer_cnv object
saveRDS(
  object = infercnv_object,
  file = paste0(out_path, time_stamp(),"infercnv_object_final.rds")
)

#######################################
# Calculate and export Infer_CNV score
#######################################


# Create a score with Extract the DE genes potentially due to CNV (count cells with more than 10% of genes altered by potential CNV)
scores = apply(infercnv_object@expr.data,2,function(x){sum(x < 0.90 | x > 1.10)/length(x)})

saveRDS(object = scores,
        file = paste0(out_path, time_stamp(),"High_CNV_scores_per_cell_10pct.rds"))


