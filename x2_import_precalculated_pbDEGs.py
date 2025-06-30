# %% =============================
import pandas as pd
import json
import os
import numpy as np


#%% ============================================
dataset_folder = "datasets/PD5D_MTG_VisiumST"

# Load the dataset metadata
data_folder = "Seurats/complex_assignment_casewilcoxauc_renamed"
files = os.listdir(data_folder)

#%% ============================================
cluster_pb_DEGs = pd.DataFrame(columns=["cluster_DE","p_val","avg_log2FC","p_val_adj","gene"])
cluster_pb_DEGs_topN = pd.DataFrame(columns=["cluster_DE","p_val","avg_log2FC","p_val_adj","gene"])
topN = 10  # Number of top DEGs to keep
for file in files:
    cell_type = file.split(".")[0]
    compare = file.split(".")[1]

    print(f"Processing {cell_type}-{compare}...")
    df = pd.read_csv(os.path.join(data_folder, file), index_col=None)
    df["cluster_DE"] = cell_type + "." + compare
    df = df.rename(columns={"pval": "p_val", "logFC": "avg_log2FC", "padj": "p_val_adj", "feature": "gene"})
    df = df[["cluster_DE", "gene", "p_val", "avg_log2FC", "p_val_adj"]]
    cluster_pb_DEGs = pd.concat([cluster_pb_DEGs, df], ignore_index=True)
    
    # get top 10 upregulated DE genes and downregulated, base on logFC
    df_up = df[df["avg_log2FC"] > 0].sort_values(by="avg_log2FC", ascending=False).head(topN)
    df_down = df[df["avg_log2FC"] < 0].sort_values(by="avg_log2FC").head(topN)
    df_topN = pd.concat([df_up, df_down], ignore_index=True)

    cluster_pb_DEGs_topN = pd.concat([cluster_pb_DEGs_topN, df_topN], ignore_index=True)

#%% ============================================
# Save the combined DataFrame to a CSV file
output_file =  dataset_folder + "/clustermarkers/cluster_pb_DEGs.csv"
cluster_pb_DEGs.to_csv(output_file, index=False)
# Save the top N DEGs DataFrame to a CSV file
output_file_topN = dataset_folder + "/clustermarkers/cluster_pb_DEGs_topN.csv"
cluster_pb_DEGs_topN.to_csv(output_file_topN, index=False)

# %% =============================================
file_ls = os.listdir(data_folder)
all_columns = []
for file in file_ls:
    cell_type = file.split(".")[0]
    compare = file.split(".")[1]

    df = pd.read_csv(os.path.join(data_folder, file), index_col=None)
    coulmn_names = df.columns.tolist()[6:]
    ## BN0737_(GLU_Neurons)THEMIS_FGF10_HC
    new_columns = [f"{col.split("_")[-2]}_{cell_type}_{col.split('_')[-1]}" for col in coulmn_names]
    all_columns.extend(new_columns)

all_columns = list(set(all_columns))

# %% ============================================
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
    file_ls = os.listdir(data_folder)
    for file in file_ls:
        cell_type = file.split(".")[0]
        compare = file.split(".")[1]

        df = pd.read_csv(os.path.join(data_folder, file), index_col=None)
        if gene in df['feature'].values:
            for col_i in df.columns[6:]:
                new_col = f"{col_i.split('_')[-2]}_{cell_type}_{col_i.split('_')[-1]}"
                gene_expr[gene][new_col] = df[df['feature'] == gene][col_i].values[0].round(2)
        else:
            continue

# %% ============================================
# Create a DataFrame from the gene expression data
gene_expr_df = pd.DataFrame(gene_expr).T
print(gene_expr_df.shape)

## assign missing values to 0
gene_expr_df = gene_expr_df.fillna(0)
print(gene_expr_df.shape)

gene_expr_df = gene_expr_df.to_csv(dataset_folder + "/clustermarkers/pb_expr_matrix_DEGs_topN.csv")
# %%
