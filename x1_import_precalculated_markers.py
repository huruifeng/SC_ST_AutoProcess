# %% =============================
import pandas as pd
import json
import os
import numpy as np


#%% ============================================
data_file = "Seurats/Jie/pseudobulk_layer_all_deseq2_contsonly_12052023.csv"

cluster_col = "layer"
gene_col = "gene"
log2fc_col = "log2FoldChange"
padj_col = "padj"
avg_expr_col = "baseMean"


dataset_folder = "datasets/PD5D_MTG_VisiumST"
meta_cluster_col = "smoothed_label_s5"
meta_sex_col = "sex"

meta_condition_col = "Condition"

output_folder = dataset_folder + "/clustermarkers"
if os.path.exists(output_folder) is False:
    os.makedirs(output_folder, exist_ok=True)



#%% ============================================
# Load the data
df = pd.read_csv(data_file, sep=",", header=0)

## fileter out the rows where the padj value is greater than 0.05
df = df[df[padj_col] <= 0.05]

## get th top 10 markers for each cluster ordered by log2fc
df = df.sort_values(by=[cluster_col, log2fc_col], ascending=[True, False])
df_top = df.groupby(cluster_col).head(10)

# only keep the relevant columns
df_top = df_top[[cluster_col, gene_col, log2fc_col, padj_col]]
# rename the columns to match the expected format
df_top = df_top.rename(columns={
    cluster_col: "cluster",
    gene_col: "gene",
    log2fc_col: "avg_log2FC",
    padj_col: "p_val_adj"
})
## values in the log2fc and padj column should be rounded to 3 decimal places
df_top["avg_log2FC"] = df_top["avg_log2FC"].round(2)
df_top["p_val_adj"] = df_top["p_val_adj"].round(2)
## save the top markers to a CSV file
df_top.to_csv(output_folder + "/cluster_markergenes_topN.csv", index=False)


#%% ============================================
# Create a dictionary to hold the markers for each cluster
# format: {cluster1: [["gene1", "log2fc1", "padj1"], ["gene2", "log2fc2", "padj2"], ...], cluster2: [["gene1", "log2fc1", "padj1"], ["gene2", "log2fc2", "padj2"], ...}
marker_genes_dict = {}
for cluster, group in df_top.groupby("cluster"):
    marker_genes_dict[cluster] = group[["gene", "avg_log2FC", "p_val_adj"]].values.tolist()

# Save the markers dictionary to a JSON file
with open(output_folder + "/cluster_markergenes_topN.json", "w") as f:
    json.dump(marker_genes_dict, f, indent=2)

#%% ============================================
pool_genes = df_top["gene"].unique().tolist()
cluster_list = df_top["cluster"].unique().tolist()
metadata = pd.read_csv(dataset_folder + "/cellspot_metadata_original.csv", index_col=0, header=0)

pct_detected = {}

marker_genes_df = pd.DataFrame()
for cluster in cluster_list:
    print("==========================")
    print("Processing cluster: ", cluster)
    pct_detected[cluster] = []
    cells_in_cluster = metadata[metadata[meta_cluster_col] == cluster]
    num_cells = len(cells_in_cluster)
    

    conditions = metadata[meta_condition_col].unique().tolist()
    sex = metadata[meta_sex_col].unique().tolist()
    for condition in conditions:
        for s in sex:
            subgroup_counts = {}
            print("Processing condition: ", condition, " and sex: ", s)
            diagnosis_sex_group = cells_in_cluster[(cells_in_cluster[meta_condition_col] == condition) & (cells_in_cluster[meta_sex_col] == s)]
            n_cells = len(diagnosis_sex_group)
            subgroup_counts["condition"] = condition
            subgroup_counts["sex"] = s
            subgroup_counts["count"] = n_cells
            pct_detected[cluster].append(subgroup_counts)  # add the subgroup counts to the list
            
    marker_genes = {}
    for gene in pool_genes:
        marker_genes[gene] = {}
        gene_expr = json.load(open(dataset_folder + "/gene_jsons/"+gene+".json"))
        gene_expr_in_cells =  list(gene_expr.keys())

        cells_with_gene_expr = [cell for cell in gene_expr_in_cells if cell in cells_in_cluster.index]
        num_cells_with_gene_expr = len(cells_with_gene_expr)

        if num_cells_with_gene_expr == 0:
            avg_expr = 0.00
        else:
            # Calculate the average expression value for the cells with gene expression and convert to float
            # avg_expr = np.mean([float(gene_expr[cell]) for cell in cells_with_gene_expr])

            # Use np.nanmean to ignore NaN values in the calculation
            # This will return NaN if all values are NaN, so we need to handle that case
            avg_expr = np.nanmean([float(gene_expr[cell]) for cell in cells_with_gene_expr])
            # If all values are NaN, set avg_expr to 0.00
            if np.isnan(avg_expr):
                avg_expr = 0.00
            else:
                avg_expr = round(avg_expr, 2)

        if avg_expr_col != "" and avg_expr_col in df.columns:
            sub_df = df[(df[cluster_col] == cluster) & (df[gene_col] == gene)]
            if sub_df.empty:
                pass  # If no sub_df, no action needed
            else:
                # Get the average expression value from the sub_df
                avg_expr = sub_df[avg_expr_col].values[0]
        else:
            pass
        marker_genes[gene]["avg_expr"] = round(avg_expr, 2)
        marker_genes[gene]["is_marker"] = gene in [g[0] for g in marker_genes_dict[cluster]] 
        marker_genes[gene]["n_expr_cells"] = num_cells_with_gene_expr
    
    marker_df = pd.DataFrame.from_dict(marker_genes, orient='index')
    marker_df[meta_cluster_col] = cluster
    marker_df["cluster_n_cells"] = num_cells

    marker_genes_df = pd.concat([marker_genes_df, marker_df], axis=0)
marker_genes_df["gene"] = marker_genes_df.index
marker_genes_df = marker_genes_df.reset_index(drop=True)
marker_genes_df = marker_genes_df[["gene", meta_cluster_col, "cluster_n_cells"] + [col for col in marker_genes_df.columns if col not in ["gene", meta_cluster_col, "cluster_n_cells"]]]  
marker_genes_df.to_csv(output_folder + "/cluster_markergenes_cellcounts.csv", index=False)

# Save the dictionary to a JSON file
with open(output_folder + "/cluster_cellcounts.json", "w") as f:
    json.dump(pct_detected, f, indent=4)    
# %%
