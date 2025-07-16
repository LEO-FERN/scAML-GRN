import os
import pandas as pd
import numpy as np

# Set the random seed for reproducibility
np.random.seed(7)

# Input directories containing the imputed.csv files for different cell types
imputed_directories = {
    "Progenitor": "~/Data/Final_Progenitor_Net",
    "Monocyte": "~/Data/Final_Monocyte_Net",
    "Dendritic": "~/Data/Final_Dendritic_Net"
}

# Directory containing the Top 1000 genes files
top_genes_dir = "~/SingleCellData"

# Output directories for filtered gene count files
filtered_output_directories = {
    "Progenitor": "~/SingleCellData/Progenitor",
    "Monocyte": "~/SingleCellData/Monocyte",
    "Dendritic": "~/SingleCellData/Dendritic"
}

# Number of samples to select for each cell type
sample_counts = {
    "Progenitor": 204,
    "Monocyte": 118,
    "Dendritic": 147
}

print("Processing...")

def filter_imputed_files(imputed_dir, top_genes_file, output_dir, cell_type, n_samples):
    """
    Filter imputed.csv files to include only the Top 1000 genes and randomly select specified number of samples.
    """
    # Read the Top 1000 genes file
    top_genes_path = os.path.join(os.path.expanduser(top_genes_dir), top_genes_file)
    try:
        top_genes = pd.read_csv(top_genes_path, sep="\t")["Gene"].tolist()
    except Exception as e:
        print(f"Error reading {top_genes_path}: {e}")
        return

    # Process all imputed.csv files in the directory
    for file in os.listdir(imputed_dir):
        if file.endswith("imputed.csv"):
            file_path = os.path.join(imputed_dir, file)

            # Load the gene count matrix
            try:
                gene_counts = pd.read_csv(file_path, header=0)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
                continue

            # Filter columns to include only Top 1000 genes
            filtered_counts = gene_counts[top_genes]

            # Randomly select the specified number of samples (rows)
            if filtered_counts.shape[0] > n_samples:
                selected_samples = np.random.choice(filtered_counts.index, size=n_samples, replace=False)
                filtered_counts = filtered_counts.loc[selected_samples]
            else:
                print(f"Warning: {file} has fewer than {n_samples} samples. Using all available samples.")

            # Write the filtered data to a new file
            output_file = os.path.join(
                os.path.expanduser(output_dir),
                file.replace("imputed.csv", "filtered.csv")
            )
            filtered_counts.to_csv(output_file, index=False)
            print(f"Filtered file saved: {output_file}")

# Process each cell type
for cell_type, imputed_dir in imputed_directories.items():
    top_genes_file = f"Top_1000_Genes_{cell_type}.txt"
    output_dir = filtered_output_directories[cell_type]
    n_samples = sample_counts[cell_type]
    filter_imputed_files(os.path.expanduser(imputed_dir), top_genes_file, os.path.expanduser(output_dir), cell_type, n_samples)
