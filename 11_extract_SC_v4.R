## This script is used to extract the data from the Seurat object.
## it will extract the data from the Seurat object and save it in a format that can be used for downstream analysis.

## Run R script as follows:
## Rscript extract_data_SC.R seurat_obj.rds output_directory

library(jsonlite)
library(Matrix)
library(data.table)
library(Seurat)
library(tidyverse)
library(jsonlite)
library(presto)

cat("===================================================\n")
# Check if the script is run with the correct number of arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Please provide the Seurat object file and output directory as arguments.")
}

# Get the arguments
seurat_obj_file <- args[1]
output_dir <- args[2]
cluster_col <- args[3]


# Check if the output directory exists, if not create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## ====================================================
# Load the Seurat object
cat("load RDS data...\n")
## Read the rds onject
seurat_obj <- readRDS(seurat_obj_file)
# Check if the Seurat object is valid
if (!inherits(seurat_obj, "Seurat")) {
  stop("The provided file is not a valid Seurat object.")
}
capture.output(str(seurat_obj), file = paste0(output_dir, "/seurat_obj_structure.txt"))

# Check if the Seurat object has the necessary assay
if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("The Seurat object does not contain the 'RNA' assay.")
}
# Check if the Seurat object has the necessary assay data
if (!"data" %in% slotNames(seurat_obj@assays$RNA)) {
  stop("The Seurat object does not contain the 'data' slot in the 'RNA' assay.")
}

# Check if the Seurat object has the necessary metadata
if (!"meta.data" %in% slotNames(seurat_obj)) {
  stop("The Seurat object does not contain the 'meta.data' slot.")
}

## Check if the Seurat object has the necessary umap embedding
if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("The Seurat object does not contain the 'umap' reduction.")
}

# Check if the cluster column exists in the metadata
if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste("The cluster column", cluster_col, "does not exist in the metadata."))
}

## ===================================================
# Extract the data
cat("Extract data...\n")

cat("Save metadata...\n")
# Save to CSV with index as the first column
metadata <- seurat_obj@meta.data
write.csv(metadata, file = paste0(output_dir, "/raw_metadata.csv"), row.names = TRUE)

## Save the metadata column names in a JSON file
cat("Save metadata feature names...\n")
# Extract metadata column names
column_names <- colnames(metadata)
# Save column names to JSON
write_json(column_names, path = file.path(output_dir, "raw_metadata_columns.json"), pretty = TRUE)

# save the umap embedding
cat("Save umap embedding...\n")
umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
colnames(umap_embeddings) <- toupper(colnames(umap_embeddings))
write.csv(umap_embeddings, file = paste0(output_dir, "/raw_umap_embeddings.csv"), row.names = TRUE)


# Extract the normalized counts
cat("Save normalized counts...\n")
## Extract normalized data
normalized_counts <- seurat_obj@assays$RNA@data  # This is a sparse matrix

# Convert sparse matrix to triplet format (long format)
long_data <- summary(normalized_counts)
# Get row (gene) and column (cell) names
long_data$Gene <- rownames(normalized_counts)[long_data$i]
long_data$Cell <- colnames(normalized_counts)[long_data$j]
long_data$Expression <- long_data$x
# Keep only necessary columns
long_data <- long_data[, c("Gene", "Cell", "Expression")]

# Filter out zero values
nonzero_data <- long_data[long_data$Expression > 0, ]

# Convert to data.table for efficiency
nonzero_data <- as.data.table(nonzero_data)
# Save to CSV with index as the first column in a fast way
fwrite(nonzero_data, file = paste0(output_dir, "/raw_normalized_counts.csv"), row.names = FALSE)


cat("Done! Data extraction is done! ^_^ ...\n")



