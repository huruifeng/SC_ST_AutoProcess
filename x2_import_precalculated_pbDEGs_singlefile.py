# %% =============================
import pandas as pd
import json
import os
import numpy as np


#%% ============================================
dataset_folder = "datasets/PD5D_MTG_VisiumST"

# Load the dataset metadata
DE_file = "Seurats/Jie/pseudobulk_de_all_deseq2_nctx_temporal_fix_12052023.csv"
pb_expr_file = "Seurats/Jie/pseudobulk_pd_layer_all_mean_lognorm_12052023.csv"

cluster_DE_col= "layer"
pval_col = "pvalue"
logFC_col = "log2FoldChange"
padj_col = "padj"
feature_col = "gene"

#%% ============================================
cluster_pb_DEGs = pd.DataFrame(columns=["cluster_DE","p_val","avg_log2FC","p_val_adj","gene"])
cluster_pb_DEGs_topN = pd.DataFrame(columns=["cluster_DE","p_val","avg_log2FC","p_val_adj","gene"])
topN = 10  # Number of top DEGs to keep

df = pd.read_csv(DE_file, index_col=0, header=0)
# Filter out rows where the p-value is greater than 0.05
df = df[df[padj_col] <= 0.05]

all_cluster_DEs = df[cluster_DE_col].unique().tolist()

for cluster_DE in all_cluster_DEs:
    print(f"Processing cluster_DE: {cluster_DE}...")
    sub_df = df[df[cluster_DE_col] == cluster_DE]
    if sub_df.empty:
        print(f"No data found for cluster_DE: {cluster_DE}. Skipping...")
        continue
    sub_df["cluster_DE"] = cluster_DE + ".CasevsControl" 
    sub_df = sub_df.rename(columns={pval_col: "p_val", logFC_col: "avg_log2FC", padj_col: "p_val_adj", feature_col: "gene"})
    sub_df = sub_df[["cluster_DE", "gene", "p_val", "avg_log2FC", "p_val_adj"]]
    cluster_pb_DEGs = pd.concat([cluster_pb_DEGs, sub_df], ignore_index=True)
    
    # get top 10 upregulated DE genes and downregulated, base on logFC
    sub_df_up = sub_df[sub_df["avg_log2FC"] > 0].sort_values(by="avg_log2FC", ascending=False).head(topN)
    sub_df_down = sub_df[sub_df["avg_log2FC"] < 0].sort_values(by="avg_log2FC").head(topN)
    sub_df_topN = pd.concat([sub_df_up, sub_df_down], ignore_index=True)

    cluster_pb_DEGs_topN = pd.concat([cluster_pb_DEGs_topN, sub_df_topN], ignore_index=True)

#%% ============================================
# Save the combined DataFrame to a CSV file
output_file =  dataset_folder + "/clustermarkers/cluster_pb_DEGs.csv"
cluster_pb_DEGs.to_csv(output_file, index=False)
# Save the top N DEGs DataFrame to a CSV file
output_file_topN = dataset_folder + "/clustermarkers/cluster_pb_DEGs_topN.csv"
cluster_pb_DEGs_topN.to_csv(output_file_topN, index=False)

# %% ============================================
pb_expr_df = pd.read_csv(pb_expr_file, index_col=0, header=0)

## colname format: SampleName_ClusterName_Condition, e.g., BN0737_TCell_HC
col_names = pb_expr_df.columns.tolist()
pb_expr_df.columns = [f"{col.split('_')[1]}_{col.split('_')[2]}_{col.split('_')[0]}" for col in col_names]


gene_pool = []
gene_expr = {}
for row in cluster_pb_DEGs_topN.itertuples():
    gene = row.gene
    print(f"Processing gene: {gene}")
    if gene not in gene_pool:
        gene_pool.append(gene)
    else:
        continue

    gene_expr[gene] = {}

    if gene in pb_expr_df.index:
        for col_i in pb_expr_df.columns:
            gene_expr[gene][col_i] = pb_expr_df.loc[gene,col_i].round(2)
    else:
        print(f"Gene {gene} not found in pb_expr_df, skipping...")

# %% ============================================
# Create a DataFrame from the gene expression data
gene_expr_df = pd.DataFrame(gene_expr).T
print(gene_expr_df.shape)

## assign missing values to 0
gene_expr_df = gene_expr_df.fillna(0)
print(gene_expr_df.shape)

gene_expr_df = gene_expr_df.to_csv(dataset_folder + "/clustermarkers/pb_expr_matrix_DEGs_topN.csv")
# %%
