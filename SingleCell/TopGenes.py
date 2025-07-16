import os
from os.path import expanduser, join
import pandas as pd

# Directories containing the adjacency matrices
directories = {
    "Progenitor": "~/BinaryFinal/Progenitor",
    "Monocyte": "~/BinaryFinal/Monocyte",
    "Dendritic": "~/BinaryFinal/Dendritic"
}

# Output directory for results
output_dir = expanduser("~/SingleCellData")
os.makedirs(output_dir, exist_ok=True)

print('Processing...')

def process_adjacency_matrices(directory, cell_type):
    """
    Process all adjacency matrix CSV files in a directory to find the top 1000 genes with the highest connectivity.
    """
    gene_connectivity = {}

    # Iterate through all CSV files in the directory
    for file in os.listdir(directory):
        if file.endswith("consensus_network.csv"):
            file_path = join(directory, file)

            # Load the adjacency matrix
            try:
                adj_matrix = pd.read_csv(file_path, header=0, index_col=None)
                adj_matrix.index = adj_matrix.columns
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
                continue

            # Calculate connectivity (sum of a single axis since the matrix is binary and symmetric)
            connectivity = adj_matrix.sum(axis=1)

            # Update global connectivity scores
            for gene, conn_value in connectivity.items():
                if gene in gene_connectivity:
                    gene_connectivity[gene] += conn_value
                else:
                    gene_connectivity[gene] = conn_value

    # Sort genes by connectivity and get the top 1000
    sorted_genes = sorted(gene_connectivity.items(), key=lambda x: x[1], reverse=True)[:1000]

    # Save results to a file
    output_file = join(output_dir, f"Top_1000_Genes_{cell_type}.txt")
    with open(output_file, "w") as f:
        f.write("Gene\tConnectivity\n")
        for gene, conn_value in sorted_genes:
            f.write(f"{gene}\t{conn_value}\n")

    print(f"Top 1000 genes with connectivity for {cell_type} saved to {output_file}")

# Process each directory
for cell_type, directory in directories.items():
    process_adjacency_matrices(expanduser(directory), cell_type)
