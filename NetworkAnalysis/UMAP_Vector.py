import os
import re
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import umap
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D

def load_and_process_network(file_path):
    df = pd.read_csv(file_path, header=0, index_col=None)
    df.index = df.columns
    return df, set(df.columns)

def get_common_genes(all_gene_sets):
    return sorted(set.intersection(*all_gene_sets))  # Sorted for consistent ordering

def extract_upper_triangle(matrix):
    return matrix.values[np.triu_indices(matrix.shape[0], k=1)]

def process_network(args):
    file_path, cell_type, common_genes = args
    network, _ = load_and_process_network(file_path)

    present_genes = [g for g in common_genes if g in network.index]
    if len(present_genes) != len(common_genes):
        missing = set(common_genes) - set(present_genes)
        raise ValueError(f"{file_path} is missing genes: {missing}")

    network = network.loc[common_genes, common_genes]
    upper_triangle = extract_upper_triangle(network)

    if np.isnan(upper_triangle).any():
        raise ValueError(f"NaNs found in upper triangle of {file_path}")

    return upper_triangle, os.path.basename(file_path), cell_type

def process_directory(directory, cell_type, common_genes):
    file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('consensus_network.csv')]
    args = [(file_path, cell_type, common_genes) for file_path in file_paths]

    with Pool(processes=10) as pool:
        results = pool.map(process_network, args)

    data, filenames, cell_types = zip(*results)
    df = pd.DataFrame(data)
    df['Filename'] = filenames
    df['Cell_Type'] = cell_types
    return df

def cluster_networks(dendritic_dir, progenitor_dir, monocyte_dir, output_plot_path, umap_coords_path):
    # Check if UMAP coordinates file already exists
    if os.path.exists(umap_coords_path):
        print(f"UMAP coordinates already saved at {umap_coords_path}. Skipping UMAP processing.")
        umap_df = pd.read_csv(umap_coords_path)
    else:
        # Find common genes
        all_gene_sets = []
        for directory in [dendritic_dir, progenitor_dir, monocyte_dir]:
            for file in os.listdir(directory):
                if file.endswith('consensus_network.csv'):
                    _, genes = load_and_process_network(os.path.join(directory, file))
                    all_gene_sets.append(genes)

        common_genes = get_common_genes(all_gene_sets)
        print(f"Number of common genes: {len(common_genes)}")

        # Process all directories
        print("Processing Dendritic cells...")
        dendritic_df = process_directory(dendritic_dir, 'Dendritic', common_genes)
        print("Processing Progenitor cells...")
        progenitor_df = process_directory(progenitor_dir, 'Progenitor', common_genes)
        print("Processing Monocyte cells...")
        monocyte_df = process_directory(monocyte_dir, 'Monocyte', common_genes)

        combined_df = pd.concat([dendritic_df, progenitor_df, monocyte_df], ignore_index=True)

        # Extract Network IDs
        combined_df['Network_ID'] = combined_df['Filename'].apply(lambda x: re.search(r'(AML\d+[A-Z]?|BM\d+)', x).group(1))

        # Run UMAP
        feature_columns = [col for col in combined_df.columns if col not in ['Filename', 'Cell_Type', 'Network_ID']]
        features = combined_df[feature_columns].values
        cell_types = combined_df['Cell_Type'].values

        print("Checking for NaNs in features:", np.isnan(features).any())
        print("Running 3D UMAP...")
        reducer = umap.UMAP(n_components=3, random_state=42)
        embedding = reducer.fit_transform(features)

        # Create DataFrame
        umap_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2', 'UMAP3'])
        umap_df['Cell_Type'] = cell_types
        umap_df['Network_ID'] = combined_df['Network_ID'].values

        # Save UMAP coordinates to file for future use
        umap_df.to_csv(umap_coords_path, index=False)
        print(f"UMAP coordinates saved to {umap_coords_path}")

    # Plot in 3D
    print("Plotting 3D UMAP...")
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection='3d')

    # Set 3D pane background color to white
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # Add faint gridlines
    ax.grid(True)
    ax.xaxis._axinfo['grid'].update(color='lightgray', linestyle='--', linewidth=0.5)
    ax.yaxis._axinfo['grid'].update(color='lightgray', linestyle='--', linewidth=0.5)
    ax.zaxis._axinfo['grid'].update(color='lightgray', linestyle='--', linewidth=0.5)

    # Create color mapping for Network_IDs
    unique_patients = sorted(umap_df['Network_ID'].unique())
    cmap = plt.cm.get_cmap('tab20', len(unique_patients))
    patient_colors = {pid: cmap(i) for i, pid in enumerate(unique_patients)}

    # Define markers for cell types
    cell_markers = {
        'Dendritic': 'o',  # Circle
        'Progenitor': 's',  # Square
        'Monocyte': '^'  # Triangle
    }

    # Plot each cell type with corresponding marker
    for cell_type, marker in cell_markers.items():
        subset = umap_df[umap_df['Cell_Type'] == cell_type]
        ax.scatter(subset['UMAP1'], subset['UMAP2'], subset['UMAP3'],
                   c=subset['Network_ID'].map(patient_colors),
                   marker=marker,
                   s=150,
                   alpha=0.7,
                   edgecolor='w',
                   linewidth=0.5,
                   facecolor='white',
                   label='_nolegend_')  # Suppress automatic legend

    # Create patient ID legend elements
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label=pid,
                   markerfacecolor=color, markersize=15)
        for pid, color in patient_colors.items()
    ]

    # Add patient legend outside plot
    ax.legend(handles=legend_elements, title='Patient IDs',
              loc='center left', bbox_to_anchor=(1.1, 0.5))

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_zlabel('UMAP3')
    ax.set_title('UMAP of Vectorized Consensus Networks (Common Genes)', pad=20, fontsize=14, fontweight='bold')

    # Remove text annotations (original code had Network_ID labels)
    #ax.grid(False)
    #fig.patch.set_facecolor('white')

    plt.savefig(output_plot_path, bbox_inches='tight')
    plt.close()

    print(f"3D UMAP plot saved to {output_plot_path}")

    #silhouette_avg = silhouette_score(embedding, cell_types)
    #print(f"The average silhouette score is: {silhouette_avg}")

def main():
    base_path = os.path.expanduser("~/BinaryFinal")
    dendritic_dir = os.path.join(base_path, "Dendritic")
    progenitor_dir = os.path.join(base_path, "Progenitor")
    monocyte_dir = os.path.join(base_path, "Monocyte")
    output_plot_path = "umap_results.png"
    umap_coords_path = "umap_coordinates.csv"  # Path to save UMAP coordinates

    cluster_networks(dendritic_dir, progenitor_dir, monocyte_dir, output_plot_path, umap_coords_path)

if __name__ == "__main__":
    main()
