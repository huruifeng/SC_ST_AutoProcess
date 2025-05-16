# %% ======================
import pandas as pd

# %% ======================
celltypemarkers = "Seurats/celltype_markers.tsv"
df = pd.read_csv(celltypemarkers, index_col=0, header=0, sep="\t")
marker_groups = df.group.unique()


# %% =====================
## read metadata
metadata = "datasets/PD5D_MTG_snRNAseq/raw_metadata.csv"
metadata_df = pd.read_csv(metadata, index_col=0, header=0)
metadata_groups = metadata_df.Complex_Assignment.unique()


# %% =====================
group_diff = set(marker_groups) - set(metadata_groups)
print("Groups in celltype_markers but not in metadata: ", group_diff)

group_diff = set(metadata_groups) - set(marker_groups)
print("Groups in metadata but not in celltype_markers: ", group_diff)



# %%
