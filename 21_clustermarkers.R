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
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 6) {
# 	stop("Please provide the correct arguments.")
# }

# # Get the arguments
# seurat_obj_file <- args[1]
# output_dir <- args[2]
# cluster_col <- args[3]
# condition_col <- args[4]
# sample_col <- args[5]
# seurat_type <- args[6]


seurat_obj_file <- "Seurats/pmdbs/pmdbs_lee_obj_updated_SN.rds"
output_dir <- "datasets/PMDBS_SN_snRNAseq"
cluster_col <- "cell_type"
condition_col <- "case"
sample_col <- "sample_id"
seurat_type <- "snrnaseq"


clustermarkers_folder = paste0(output_dir, "/clustermarkers")
if (!dir.exists(clustermarkers_folder)) {
  dir.create(clustermarkers_folder, recursive = TRUE)
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
if(seurat_type == "scrnaseq" | seurat_type == "snrnaseq"){
	  # Check if the Seurat object has the necessary assay
	if (!"RNA" %in% names(seurat_obj@assays)) {
		stop("The Seurat object does not contain the 'RNA' assay.")
	}
	# Check if the Seurat object has the necessary assay data
	if (!"data" %in% names(seurat_obj@assays$RNA@layers)) {
		stop("The Seurat object does not contain the 'data' slot in the 'RNA' assay.")
	}
} else if (seurat_type == "snatacseq" | seurat_type == "scatacseq") {
	# Check if the Seurat object has the necessary assay
	if (!"ATAC" %in% names(seurat_obj@assays)) {
		stop("The Seurat object does not contain the 'ATAC' assay.")
	}
	# Check if the Seurat object has the necessary assay data
	if (!"counts" %in% slotNames(seurat_obj@assays$ATAC)) {
		stop("The Seurat object does not contain the 'counts' slot in the 'ATAC' assay.")
	}
} else if (seurat_type == "visiumst") {
	# Check if the Seurat object has the necessary assay
	if (!"Spatial" %in% names(seurat_obj@assays)) {
		stop("The Seurat object does not contain the 'Spatial' assay.")
	}
	# Check if the Seurat object has the necessary assay data
	if (!"data" %in% slotNames(seurat_obj@assays$Spatial)) {
		stop("The Seurat object does not contain the 'data' slot in the 'Spatial' assay.")
	}
}


# Check if the Seurat object has the necessary metadata
if (!"meta.data" %in% slotNames(seurat_obj)) {
  stop("The Seurat object does not contain the 'meta.data' slot.")
}

# Check if the cluster column exists in the metadata
if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste("The cluster column", cluster_col, "does not exist in the metadata."))
}


## ===========================================================
# cell type specific markers
print("=====================================================")
print("Extracting cell type specific markers...")

# Extract cell type specific markers
cell_type_markers <- FindAllMarkers(seurat_obj, group.by =  cluster_col)
# Convert to data.table
cell_type_markers_dt <- as.data.table(cell_type_markers)
# Save to CSV
fwrite(cell_type_markers_dt, paste0(clustermarkers_folder,"/cluster_FindAllMarkers.csv"), row.names = FALSE)

## ============================================================
# Ecalcaulate differential expression within each cell type between the conditions # nolint
print("=====================================================")
print("Calculating differential expression within each cell type...")

# Define the cell types
cell_types <- unique(seurat_obj@meta.data[[cluster_col]])
# Initialize an empty list to store results
de_results_list <- list()
de_results_topN_list <- list()

# Loop through each cell type
for (cell_type in cell_types) {
	# Print the current cell type # nolint: whitespace_linter, indentation_linter.
	print(paste("Processing cell type:", cell_type))
	# Subset the Seurat object to the current cell type # nolint
	subset_obj <- seurat_obj[, seurat_obj@meta.data[[cluster_col]] == cell_type]

	# Calculate differential expression between conditions
	condition_ls <- unique(seurat_obj@meta.data[[condition_col]])
	# if there are more than 2 conditions, compare all possible combinations
	if (length(condition_ls) > 2) {
		combinations <- combn(condition_ls, 2)
		for (i in 1:ncol(combinations)) {
			c1 <- combinations[1, i]
			c2 <- combinations[2, i]
			de_results <- FindMarkers(subset_obj, ident.1 = c1, ident.2 = c2, group.by = condition_col)
			## add a column for the gene names
			de_results$gene <- rownames(de_results)
			
			## filter out genes with padj > 0.05
			# de_results <- de_results[de_results$p_val_adj < 0.05, ]
			
			## get top 10 upregulated DE genes and downregulated, base on logFC
			de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
		
			# Store the results in the list
			de_results_list[[paste(cell_type, paste(c1, c2, sep = "vs"), sep = ".")]] <- de_results
			de_results_topN_list[[paste(cell_type, paste(c1, c2, sep = "vs"), sep = ".")]] <- de_results_topN
		}
	} else {
		de_results <- FindMarkers(subset_obj, ident.1 = condition_ls[1], ident.2 = condition_ls[2], group.by = condition_col)
		## add a column for the gene names
		de_results$gene <- rownames(de_results)

		## filter out genes with padj > 0.05
		# de_results <- de_results[de_results$p_val_adj < 0.05, ]
		
		## get top 10 upregulated DE genes and downregulated, base on logFC
		de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
		
		de_results_list[[paste(cell_type, paste(condition_ls[1], condition_ls[2], sep = "vs"), sep = ".")]] <- de_results
		de_results_topN_list[[paste(cell_type, paste(condition_ls[1], condition_ls[2], sep = "vs"), sep = ".")]] <- de_results_topN

	}
}
# Convert the list to a data.table
de_results_dt <- rbindlist(de_results_list, idcol = "cluster_DE")
de_results_topN_dt <- rbindlist(de_results_topN_list, idcol = "cluster_DE")
# Save to CSV
fwrite(de_results_dt, paste0(clustermarkers_folder, "/cluster_DEGs.csv"), row.names = FALSE)
fwrite(de_results_topN_dt, paste0(clustermarkers_folder, "/cluster_DEGs_topN.csv"), row.names = FALSE)

## ============================================================
print("=====================================================")
# pseudo-bulk DE analysis in each cell type
print("Calculating pseudo-bulk differential expression within each cell type...")
# Create a new Seurat object for pseudo-bulk analysis
# pb_obj <- AggregateExpression(seurat_obj, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("sample_id", "MajorCellTypes", "case"))
if(seurat_type == "scrnaseq" | seurat_type == "snrnaseq"){
	## use the RNA assay
	assay = "RNA"
} else if (seurat_type == "snATACseq") {
	## use the ATAC assay
	assay = "ATAC"
} else if (seurat_type == "visiumst") {
	## use the Spatial assay
	assay = "Spatial"
} else {
	stop("Invalid seurat_type. Please choose from 'scrnaseq', 'snrnaseq', 'snATACseq', or 'visiumst'.")
}
## use the counts assay
pb_obj <- PseudobulkExpression(seurat_obj,assays = assay, layer = "counts", method= "aggregate", return.seurat = T, group.by = c(sample_col, cluster_col, condition_col))

pb_obj@meta.data[[cluster_col]] <- gsub("-", "_", pb_obj@meta.data[[cluster_col]])
pb_obj@meta.data[["orig.ident"]] <- gsub("-", "_", pb_obj@meta.data[["orig.ident"]])

metadata = pb_obj@meta.data
## rename sample_id to sampleId
colnames(metadata)[colnames(metadata) == sample_col] <- "sampleId"
colnames(metadata)[colnames(metadata) == condition_col] <- "condition"

write.csv(metadata, paste0(clustermarkers_folder, "/metadata_sample_cluster_condition.csv"), row.names = FALSE)

expr_matrix <- GetAssayData(pb_obj, assay = assay, layer = "data")
colnames(expr_matrix) <- gsub("-", "_", colnames(expr_matrix))
write.csv(expr_matrix, paste0(clustermarkers_folder, "/pb_expr_matrix.csv"), row.names = TRUE)


# Define the cell types
cell_types <- unique(pb_obj@meta.data[[cluster_col]])

# Initialize an empty list to store results
pseudo_bulk_list <- list()
pseudo_bulk_topN_list <- list()
# Loop through each cell type
for (cell_type in cell_types) {
	# Print the current cell type # nolint: whitespace_linter, indentation_linter.
	print(paste("Processing cell type:", cell_type))
	# Subset the Seurat object to the current cell type # nolint
	subset_obj <- pb_obj[, pb_obj@meta.data[[cluster_col]] == cell_type]

	# Calculate differential expression between conditions
	# Here, we assume that the conditions are stored in the "Condition" metadata column 
	# You may need to adjust this based on your actual metadata structure
	condition_ls <- unique(subset_obj@meta.data[[condition_col]])
	# if there are more than 2 conditions, compare all possible combinations
	if (length(condition_ls) > 2) {
		combinations <- combn(condition_ls, 2)
		for (i in 1:ncol(combinations)) {
			c1 <- combinations[1, i]
			c2 <- combinations[2, i]
			de_results <- FindMarkers(subset_obj, ident.1 = c1, ident.2 = c2, group.by = condition_col, test.use = "wilcox")
			## add a column for the gene names
			de_results$gene <- rownames(de_results)
			
			## filter out genes with padj > 0.05
			# de_results <- de_results[de_results$p_val_adj < 0.1, ]

			## get top 10 upregulated DE genes and downregulated, base on logFC
			de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
			
			# Store the results in the list
			pseudo_bulk_list[[paste(cell_type, paste(c1, c2, sep = "vs"), sep = ".")]] <- de_results
			pseudo_bulk_topN_list[[paste(cell_type, paste(c1, c2, sep = "vs"), sep = ".")]] <- de_results_topN
		}
	} else {
		de_results <- FindMarkers(subset_obj, ident.1 = condition_ls[1], ident.2 = condition_ls[2], group.by = condition_col,test.use = "wilcox")
		## add a column for the gene names
		de_results$gene <- rownames(de_results)
		
		## filter out genes with padj > 0.05
		# de_results <- de_results[de_results$p_val_adj < 0.05, ]
		
		## get top 10 upregulated DE genes and downregulated, base on logFC
		de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
		
		# Store the results in the list
		pseudo_bulk_list[[paste(cell_type, paste(condition_ls[1], condition_ls[2], sep = "vs"), sep = ".")]] <- de_results
		pseudo_bulk_topN_list[[paste(cell_type, paste(condition_ls[1], condition_ls[2], sep = "vs"), sep = ".")]] <- de_results_topN
	}
}
# Convert the list to a data.table
pseudo_bulk_dt <- rbindlist(pseudo_bulk_list, idcol = "cluster_DE")
pseudo_bulk_topN_dt <- rbindlist(pseudo_bulk_topN_list, idcol = "cluster_DE")
# Save to CSV
fwrite(pseudo_bulk_dt, paste0(clustermarkers_folder, "/cluster_pb_DEGs.csv"), row.names = FALSE)
fwrite(pseudo_bulk_topN_dt, paste0(clustermarkers_folder, "/cluster_pb_DEGs_topN.csv"), row.names = FALSE)

pooled_topN_DEGs = pseudo_bulk_topN_dt$gene
## remove duplicates
pooled_topN_DEGs <- unique(pooled_topN_DEGs)
## subset expr_matrix to only include pooled_topN_DEGs
expr_matrix_pooled_topN_DGEs <- expr_matrix[pooled_topN_DEGs, ]
## save the pooled_topN_DEGs expression matrix
write.csv(expr_matrix_pooled_topN_DGEs, paste0(clustermarkers_folder, "/pb_expr_matrix_DEGs_topN.csv"), row.names = TRUE)

cat("Done! cluster markers is done! ^_^ ...\n")



