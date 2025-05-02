library(jsonlite)
library(Matrix)
library(data.table)
library(Seurat)
library(presto)
library(tidyverse)

print("load RDS data...")
## Read the rds onject
seurat_obj_file = "ST_data/comb_seurat_smoothed_V4_cell2loc_01222025.rds"
seurat_obj <- readRDS(seurat_obj_file)

sample_ls = c("BN1827", "BN1839", "BN1424", "BN1535", "BN1726", "BN1822", "BN0662", "BN1817", "BN1762", "BN1076")

subset_obj <- subset(seurat_obj, subset = (sample_id %in% sample_ls))

saveRDS(subset_obj, file = "SC_data/FinalOSR_6HC_MTG_10samples.rds")