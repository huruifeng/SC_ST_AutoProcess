import toml
import os
import subprocess

# Check if the toml file exists
toml_file = "sc_dataset_info.toml"
if not os.path.exists(toml_file):
    print(f"File {toml_file} not found.")
    exit(1)

with open(toml_file, "r") as f:
    toml_data = toml.load(f)
    # print(toml_data)

# Check if the toml file has the required keys
seurat = toml_data.get("seurat", {})
if not seurat:
    print("No seurat info found in the toml file.")
    exit(1)

seurat_file = seurat.get("seurat_file", "")
if not seurat_file:
    print("No seurat file found in the toml file.")
    exit(1)

seurat_path = f"Seurat/{seurat_file}"
if not os.path.exists(seurat_path):
    print(f"Seurat file {seurat_path} not found.")
    exit(1)

dataset_type = seurat.get("datatype", "")
if not dataset_type:
    print("No dataset type found in the toml file.")
    exit(1)

if dataset_type.lower() not in ["scrnaseq", "visiumst"]:
    print(f"Invalid dataset type {dataset_type} found in the toml file.")
    exit(1)

## Check if the toml file has the dataset keys
dataset_info = toml_data.get("dataset_info", {})
if not dataset_info:
    print("No dataset info found in the toml file.")
    exit(1)

dataset_name = dataset_info.get("dataset_name", "")
if not dataset_name:
    print("No dataset name found in the toml file.")
    exit(1)

selected_features = dataset_info.get("selected_features", "")
if not selected_features:
    print("No selected features found in the toml file.")
    exit(1)
sample_id_column = dataset_info.get("sample_id_column", "")
if not sample_id_column:
    print("No sample id column found in the toml file.")
    exit(1)
major_cluster_column = dataset_info.get("major_cluster_column", "")
if not major_cluster_column:
    print("No major cluster column found in the toml file.")
    exit(1)
condition_column = dataset_info.get("condition_column", "")
if not condition_column:
    print("No condition column found in the toml file.")
    exit(1)

## run the R script
dataset_path = f"datasets/{dataset_name}"
 ## check if dataset exists
if not os.path.exists(dataset_path):
    os.makedirs(dataset_path)

log_file = open(f"{dataset_path}/extract_seurat_output.log", "w")
if dataset_type.lower() in ["scrnaseq", "snrnaseq"]:
    subprocess.run(
        ["Rscript", "extract_SC.R", f"Seurats/{seurat}", dataset_path],
        stdout=log_file,
        stderr=log_file,
    )
elif dataset_type.lower() in ["visiumst"]:
    subprocess.run(
        ["Rscript", "extract_Visium.R", f"Seurats/{seurat}", dataset_path],
        stdout=log_file,
        stderr=log_file,
    )
else:
    print(f"Invalid dataset type {dataset_type} found in the toml file.")
    exit(1)
log_file.close()


log_file = open(f"{dataset_path}/prepare_meta_output.log", "w")
if dataset_info["seurat"]["datatype"].lower() in ["scrnaseq", "snrnaseq"]:
    subprocess.Popen(
        ["python3", "ename_meta_SC.py",dataset_path, ",".join(selected_features), sample_id_column, major_cluster_column, condition_column],
        stdout=log_file,
        stderr=log_file,
    )
elif dataset_info["seurat"]["datatype"].lower() in ["visiumst"]:
    subprocess.Popen(
        ["python3", "rename_meta_Visium.py",dataset_path, ",".join(selected_features), sample_id_column, major_cluster_column, condition_column],
        stdout=log_file,
        stderr=log_file,
    )
else:
    print(f"Invalid dataset type {dataset_type} found in the toml file.")
    exit(1)
log_file.close()



