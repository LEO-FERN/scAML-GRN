import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import math

def create_heatmaps(input_file, output_dir, max_pathways_per_plot=30):
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Get unique cell types
    cell_types = df['CellType'].unique()

    # Determine if special chunking is needed (for C5 file)
    is_c5 = 'C5' in Path(input_file).stem

    for cell_type in cell_types:
        cell_data = df[df['CellType'] == cell_type]
        heatmap_data = cell_data.set_index('PatientID').drop(columns=['CellType'])
        heatmap_data = heatmap_data.T
        heatmap_data = heatmap_data.loc[~(heatmap_data == 0).all(axis=1)]

        if is_c5:
            # Split into chunks of max_pathways_per_plot
            total_pathways = heatmap_data.shape[0]
            num_chunks = math.ceil(total_pathways / max_pathways_per_plot)

            for i in range(num_chunks):
                chunk_data = heatmap_data.iloc[i * max_pathways_per_plot:(i + 1) * max_pathways_per_plot]

                fig_width = max(10, len(chunk_data.columns) * 0.5)
                fig_height = max(8, len(chunk_data.index) * 0.4)
                plt.figure(figsize=(fig_width, fig_height))

                sns.heatmap(chunk_data, cmap=['white', 'red'], vmin=0, vmax=1,
                            linewidths=0.5, linecolor='gray', cbar=False)  # Legend removed

                #plt.title(f'Pathway Activation - {cell_type} (Chunk {i + 1})', pad=20)
                plt.xlabel('Patient ID')
                plt.ylabel('Pathways')
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout(rect=[0.2, 0, 1, 1])

                filename = f"{Path(input_file).stem}_{cell_type.replace(' ', '_')}_chunk{i + 1}_heatmap.svg"
                output_path = os.path.join(output_dir, filename)
                plt.savefig(output_path, format='svg', dpi=300, bbox_inches='tight')
                plt.close()
        else:
            fig_width = max(10, len(heatmap_data.columns) * 0.5)
            fig_height = max(8, len(heatmap_data.index) * 0.4)
            plt.figure(figsize=(fig_width, fig_height))

            sns.heatmap(heatmap_data, cmap=['white', 'red'], vmin=0, vmax=1,
                        linewidths=0.5, linecolor='gray', cbar=False)  # Legend removed

            #plt.title(f'Pathway Activation - {cell_type}', pad=20)
            plt.xlabel('Patient ID')
            plt.ylabel('Pathways')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout(rect=[0.2, 0, 1, 1])

            filename = f"{Path(input_file).stem}_{cell_type.replace(' ', '_')}_heatmap.svg"
            output_path = os.path.join(output_dir, filename)
            plt.savefig(output_path, format='svg', dpi=300, bbox_inches='tight')
            plt.close()

    print(f"Processed {input_file} - created heatmaps for {len(cell_types)} cell types" +
          (f", each in {num_chunks} chunks" if is_c5 else ""))

def main():
    input_dir = os.path.expanduser('~/EnrichmentData')
    output_dir = os.path.expanduser('~/ThesisCode/EnrichmentAnalysis')
    os.makedirs(output_dir, exist_ok=True)

    input_files = [
        'pathway_matrix_C2.csv',
        'pathway_matrix_C5.csv',
        'pathway_matrix_C6.csv',
        'pathway_matrix_H.csv'
    ]

    for file in input_files:
        input_path = os.path.join(input_dir, file)
        if os.path.exists(input_path):
            create_heatmaps(input_path, output_dir)
        else:
            print(f"File not found: {input_path}")

if __name__ == "__main__":
    main()
