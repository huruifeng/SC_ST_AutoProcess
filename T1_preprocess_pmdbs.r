library(Seurat)

lee <- readRDS("Seurats/pmdbs_sc_rnaseq_cohort_analysis_team-lee.final_HC_all_genes.rds")
lee$brain_region[lee$brain_region == "Substantia_Nigra "]<- "Substantia_Nigra"
lee <- NormalizeData(lee, normalization.method = "LogNormalize", scale.factor = 10000)
lee[["umap"]] <- lee[["umap.original"]]
lee[["umap.original"]] <- NULL

## add a column "case" to the metadata, based on the primay diagnosis
lee$case <- ifelse(lee$primary_diagnosis == "Idiopathic PD", "PD", "Control")
# lee$case <- dplyr::case_when(
#   lee$primary_diagnosis == "Idiopathic PD" ~ "PD",
#   lee$primary_diagnosis == "No PD nor other neurological disorder" ~ "Control",
#   TRUE ~ "Other"
# )

lee_PMDBS<- read.csv("Seurats/pmdbs/metadata_release_PMDBS.csv", row.names = 1)
lee$sample_id = lee_PMDBS$sample_id[match(lee$ASAP_sample_id, lee_PMDBS$ASAP_sample_id)] 

capture.output(str(lee), file = paste0("Seurats/pmdbs/", "lee_obj_structure.txt"))

saveRDS(lee, "Seurats/pmdbs_lee_obj_updated.rds")






# ============================================
metadata<- read.csv("Seurats/pmdbs/ASAP_team_Lee/pmdbs_sc_rnaseq_cohort_analysis_team-lee.final_metadata.csv")
lee_subject_metadata<- read.csv("Seurats/pmdbs/ASAP_team_Lee/metadata_release_SUBJECT.csv", row.names = 1)
lee_PMDBS<- read.csv("Seurats/pmdbs/ASAP_team_Lee/metadata_release_PMDBS.csv", row.names = 1)

cells.meta<- lee@meta.data
cells.meta<- cells.meta%>%
  tibble::rownames_to_column("barcode")

cells.meta$subject_id <- sub("(_[^_]+){2}$", "", cells.meta$sample)
# add diagnosis info from subject metadata csv:
cells.meta<- merge(cells.meta, lee_subject_metadata[, !(names(lee_subject_metadata) %in% "subject_id")],
                   by.x = "subject_id", by.y = "ASAP_subject_id",
                   all.x = TRUE)
# add brain region from PMDBS metadata csv:
cells.meta$ASAP_sample_id<- sub("_[^_]+$", "", cells.meta$sample) # create this column to match with the PMDBS sheet

cells.meta<- merge(cells.meta, lee_PMDBS[, !(names(lee_PMDBS) %in% c("sample_id", "ASAP_team_id", "ASAP_dataset_id"))],
                  by = "ASAP_sample_id",
                   all.x = TRUE)

cells.meta<- cells.meta%>%
  tibble::column_to_rownames("barcode")

# check rowname/barcode order!!
if(!all(rownames(lee@meta.data) == rownames(cells.meta))){
  cells.meta<- cells.meta[match(rownames(lee@meta.data), rownames(cells.meta)),]
}

lee<- AddMetaData(lee, metadata = cells.meta)
