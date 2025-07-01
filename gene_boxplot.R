library(jsonlite)
library(Matrix)
library(data.table)
library(Seurat)
library(presto)
library(tidyverse)
library(ggplot2)

print("load RDS data...")
## Read the rds onject
seurat_obj_file = "Seurats/Jacobs/DietFullIntegration_Oct2023_FinalOSR_6HC_MTG_FebCA_MajorMarkersUpdated2.rds"
seurat_obj <- readRDS(seurat_obj_file)

## Option 1:
VlnPlot(seurat_obj, features = "LITAF", group.by = "MajorCellTypes", pt.size = 0) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Option 2:
# Fetch expression and metadata
df <- FetchData(seurat_obj, vars = c("LITAF", "Complex_Assignment", "MajorCellTypes"))
df <- df[df$LITAF > 0, ]

ggplot(df, aes(x = Complex_Assignment, y = LITAF)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Sub Cell Type") +
  ylab("LITAF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save with custom width and height (in inches)
ggsave("boxplot_LITAF_without_0s_subcelltype.pdf", width = 16, height = 5, units = "in", dpi = 300)


## Option 3: Multiple Genes
genes <- c("LITAF", "MS4A1")
df <- FetchData(seurat_obj, vars = c(genes, "MajorCellTypes")) %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  filter(expression > 0)  # remove 0s

ggplot(df, aes(x = MajorCellTypes, y = expression)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
