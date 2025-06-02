# %% ==============================
import pandas as pd
import json
from collections import defaultdict
import os
import numpy as np


# %% ============================================================================
## top 10 marker genes in each cell type
dataset_folder = "datasets/PD5D_MTG_snRNAseq"
# dataset_folder = "datasets/PD5D_MTG_VisiumST"

with open(f"{dataset_folder}/clustermarkers/cluster_markergenes.json", "r") as f:
    cluster_markergenes = json.load(f)
cluster_mapping = {(")").join(s.split(')')[1:]): s for s in cluster_markergenes}


# %% ============================================================================
data_file = "Seurats/PrestoFindAllMarkersTop.tsv"
df = pd.read_csv(data_file, sep="\t", header=0)

# check if the value in the cluster column are all in the cluster_mapping
if not all(df["group"].isin(cluster_mapping.keys())):
    print("Some values in the cluster column are not in the cluster_mapping.")
    print("These values will be kept as is.")
    print("Values not in the cluster_mapping:", df[~df["group"].isin(cluster_mapping.keys())]["group"].unique())

# %% ============================================================================
## if the vvalue in the cluster column is not in the cluster_mapping, keep it as is
df["group"] = df["group"].apply(lambda x: cluster_mapping[x] if x in cluster_mapping else x)


# %%
