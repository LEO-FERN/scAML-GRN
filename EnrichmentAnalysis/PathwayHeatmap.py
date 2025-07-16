import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.colors import ListedColormap

#Either '1000' or 'important' to select dataset
def plot_combined_heatmaps(input_dir, output_dir, prefix): 
    # Expand home directory shortcuts
    input_dir = os.path.expanduser(input_dir)
    output_dir = os.path.expanduser(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Define cluster processing order and custom colormap
    cluster_order = ['H', 'C2', 'C5', 'C6']
    cmap = ['white', 'red']

    # Find and sort relevant files
    files = glob.glob(os.path.join(input_dir, f"{prefix}_pathway*.csv"))
    files = sorted(files, key=lambda x: cluster_order.index(x.split('_')[-1].split('.')[0]))

    # Process each cluster file
    combined_df = pd.DataFrame()
    for file in files:
        cluster = file.split('_')[-1].split('.')[0]
        df = pd.read_csv(file).set_index('CellType')
        df.columns = [f"{cluster}_{col}" for col in df.columns]  # Add cluster ID to pathways
        combined_df = pd.concat([combined_df, df], axis=1)

    # Determine figure size dynamically based on data dimensions
    fig_width = max(10, len(combined_df.columns) * 0.4)
    fig_height = max(8, len(combined_df.index) * 0.5)

    # Create unified heatmap
    plt.figure(figsize=(fig_width, fig_height))
    sns.heatmap(
        combined_df.T,  # Transpose for pathways as rows
        cmap=cmap,
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor='gray',
        cbar=False
    )

    # Formatting
    #plt.title(f'Combined Pathway Activation - {prefix.capitalize()} Set', pad=20)
    plt.xlabel('Cell Types')
    plt.ylabel('Pathways')
    plt.xticks(rotation=45, ha='right')

    # Save output as high-quality SVG
    output_path = os.path.join(output_dir, f"{prefix}_combined_heatmap.svg")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, format='svg')
    plt.close()


# Example usage:
input_dir = "~/EnrichmentData"
output_dir = "~/ThesisCode/EnrichmentAnalysis"

# Generate combined heatmaps for both datasets
plot_combined_heatmaps(input_dir, output_dir, '1000')
plot_combined_heatmaps(input_dir, output_dir, 'important')
