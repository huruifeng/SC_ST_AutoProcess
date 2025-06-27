import toml
import os
import json
import subprocess

## ==============================================================
print("Loading toml file...")
# Check if the toml file exists
toml_file = "00_sc_dataset_info.toml"
if not os.path.exists(toml_file):
    print(f"File {toml_file} not found.")
    exit(1)

with open(toml_file, "r") as f:
    toml_data = toml.load(f)
    # print(toml_data)

## ==============================================================
print("Checking toml file...[Seurat]")
# Check if the toml file has the required keys
seurat = toml_data.get("seurat", {})
if not seurat:
    print("No seurat info found in the toml file.")
    exit(1)

seurat_file = seurat.get("seurat_file", "")
if not seurat_file:
    print("No seurat file found in the toml file.")
    exit(1)

seurat_path = f"Seurats/{seurat_file}"
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

## ==============================================================
print("Checking toml file...[Dataset]")
## Check if the toml file has the dataset keys
dataset = toml_data.get("dataset", {})
if not dataset:
    print("No dataset info found in the toml file.")
    exit(1)

dataset_name = dataset.get("dataset_name", "")
if not dataset_name:
    print("No dataset name found in the toml file.")
    exit(1)

## ==============================================================
print("Checking toml file...[Meta Features]")
meta_features = toml_data.get("meta_features", {})
if not meta_features:
    print("No meta features found in the toml file.")
    exit(1)

selected_features = meta_features.get("selected_features", [])
if not selected_features:
    print("No selected features found in the toml file.")
    exit(1)

sample_id_column = meta_features.get("sample_id_column", "")
if not sample_id_column:
    print("No sample id column found in the toml file.")
    exit(1)
main_cluster_column = meta_features.get("main_cluster_column", "")
if not main_cluster_column:
    print("No main cluster column found in the toml file.")
    exit(1)
condition_column = meta_features.get("condition_column", "")
if not condition_column:
    print("No condition column found in the toml file.")
    exit(1)

## ==============================================================
## save a copy of the toml file in the dataset folder
dataset_path = f"datasets/{dataset_name}"
 ## check if dataset exists
if not os.path.exists(dataset_path):
    os.makedirs(dataset_path)

## copy the toml file to the dataset folder
with open(f"{dataset_path}/dataset_info.toml", "w") as f:
    toml.dump(toml_data, f)

## ==============================================================
print("==================================================")
print("Running R script...Extract Seurat data...")
## run the R script
with open(f"{dataset_path}/extract_seurat_output.log", "w") as log_file:
    if dataset_type.lower() in ["scrnaseq", "snrnaseq"]:
        process = subprocess.Popen(
            ["Rscript", "11_extract_SC.R", seurat_path, dataset_path, main_cluster_column],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # get string output, not bytes
        )
    elif dataset_type.lower() in ["visiumst"]:
        process = subprocess.Popen(
            ["Rscript", "11_extract_Visium.R", seurat_path, dataset_path, main_cluster_column],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # get string output, not bytes
        )
    else:
        print(f"Invalid dataset type {dataset_type} found in the toml file.")
        exit(1)

    # Stream output live
    for line in process.stdout:
        print(line, end="")         # print to terminal
        log_file.write(line)        # write to log file
        log_file.flush()            # ensure it’s written immediately

    process.wait()


## ==============================================================
print("==================================================")
print("Running R script...Cell type markers...")
## run the R script
with open(f"{dataset_path}/cluster_markers.log", "w") as log_file:
    if dataset_type.lower() in ["scrnaseq", "snrnaseq","visiumst"]:
        process = subprocess.Popen(
            ["Rscript", "21_clustermarkers.R", seurat_path, dataset_path, main_cluster_column, condition_column, sample_id_column, dataset_type.lower()],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # get string output, not bytes
        )
    else:
        print(f"Invalid dataset type {dataset_type} found in the toml file.")
        exit(1)

    # Stream output live
    for line in process.stdout:
        print(line, end="")         # print to terminal
        log_file.write(line)        # write to log file
        log_file.flush()            # ensure it’s written immediately

    process.wait()

## ==============================================================
print("=================================================")
print("Running python script...Prepare meta data...")
## run the python script
print("Selected features:", ",".join(selected_features))
with open(f"{dataset_path}/prepare_meta_output.log", "w") as log_file:
    if dataset_type.lower() in ["scrnaseq", "snrnaseq"]:
        process_py = subprocess.Popen(
            ["python3", "31_rename_meta_SC.py",dataset_path, ",".join(selected_features), sample_id_column, main_cluster_column, condition_column],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # get string output, not bytes
        )
    elif dataset_type.lower() in ["visiumst"]:
        process_py = subprocess.Popen(
            ["python3", "31_rename_meta_Visium.py",dataset_path, ",".join(selected_features), sample_id_column, main_cluster_column, condition_column],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # get string output, not bytes
        )
    else:
        print(f"Invalid dataset type {dataset_type} found in the toml file.")
        exit(1)
    # Stream output live
    for line in process_py.stdout:
        print(line, end="")         # print to terminal
        log_file.write(line)        # write to log file
        log_file.flush()            # ensure it’s written immediately

    process_py.wait()

## ==============================================================
print("=================================================")
print("Running python script...Prepare cell count data...")
## run the python script
with open(f"{dataset_path}/prepare_cellcount_output.log", "w") as log_file:
    if dataset_type.lower() in ["scrnaseq", "snrnaseq","visiumst"]:

        process_py = subprocess.Popen(
            ["python3", "41_clustermarkers_postprocess.py", dataset_path, main_cluster_column, "sex"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # get string output, not bytes
        )

    else:
        print(f"Invalid dataset type {dataset_type} found in the toml file.")
        exit(1)

    # Stream output live
    for line in process_py.stdout:
        print(line, end="")         # print to terminal
        log_file.write(line)        # write to log file
        log_file.flush()            # ensure it’s written immediately

    process_py.wait()
## ==============================================================
print("Done! All scripts executed successfully.")




