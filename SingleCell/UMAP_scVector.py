import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import umap
from sklearn.preprocessing import StandardScaler
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D


# Helper function definitions
def process_network(file_path):
    adj_matrix = pd.read_csv(file_path, header=0, index_col=None).values
    upper_tri = adj_matrix[np.triu_indices(len(adj_matrix), k=1)]
    return upper_tri


def process_patient(args):
    patient_dir, cell_type = args
    networks = []
    for file in os.listdir(patient_dir):
        if file.endswith('.csv'):
            file_path = os.path.join(patient_dir, file)
            networks.append(process_network(file_path))
    return np.array(networks), os.path.basename(patient_dir), cell_type


def process_cell_type(cell_type_dir):
    patient_dirs = [os.path.join(cell_type_dir, d) for d in os.listdir(cell_type_dir)
                    if os.path.isdir(os.path.join(cell_type_dir, d))]
    cell_type = os.path.basename(cell_type_dir)
    args = [(patient_dir, cell_type) for patient_dir in patient_dirs]

    with Pool(processes=20) as pool:
        results = pool.map(process_patient, args)

    data, patient_ids, cell_types = zip(*results)
    flat_data = np.vstack(data)
    flat_patient_ids = np.repeat(patient_ids, [len(d) for d in data])

    df = pd.DataFrame(flat_data)
    df['Patient_ID'] = flat_patient_ids
    df['Cell_Type'] = cell_types[0]
    return df


def load_data(cell_type_dir):
    """Loads data and returns feature matrix (X) and patient labels (y)."""
    df = process_cell_type(cell_type_dir)
    X = df.drop(['Patient_ID', 'Cell_Type'], axis=1).values
    y = df['Patient_ID'].values
    return X, y


# New functions for consistent coloring
def get_all_patient_ids(base_dir, cell_types):
    """Collect all unique patient IDs across all cell types."""
    all_patient_ids = set()
    for cell_type in cell_types:
        cell_type_dir = os.path.join(base_dir, cell_type)
        patient_dirs = [d for d in os.listdir(cell_type_dir)
                        if os.path.isdir(os.path.join(cell_type_dir, d))]
        all_patient_ids.update(patient_dirs)
    return sorted(all_patient_ids)


def create_global_color_mapping(patient_ids):
    """Generate consistent color palette for all patients using Matplotlib's tab20."""
    cmap = plt.cm.get_cmap('tab20', len(patient_ids))
    return {pid: cmap(i) for i, pid in enumerate(patient_ids)}


# Modified plotting function
def run_umap_3d(cell_type_dir, cell_type, global_color_map, n_neighbors=15, min_dist=0.1):
    """Create 3D UMAP plot with consistent colors and formatting."""
    X, y = load_data(cell_type_dir)

    # Data processing
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    umap_model = umap.UMAP(n_components=3, n_neighbors=n_neighbors,
                           min_dist=min_dist, random_state=42)
    X_umap = umap_model.fit_transform(X_scaled)

    # Plot setup
    fig = plt.figure(figsize=(12, 8), facecolor='white')
    ax = fig.add_subplot(111, projection='3d', facecolor='white')

    # Set 3D pane background color to white
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # Add faint gridlines
    ax.grid(True)
    ax.xaxis._axinfo['grid'].update(color='lightgray', linestyle='--', linewidth=0.5)
    ax.yaxis._axinfo['grid'].update(color='lightgray', linestyle='--', linewidth=0.5)
    ax.zaxis._axinfo['grid'].update(color='lightgray', linestyle='--', linewidth=0.5)

    ax.view_init(elev=20, azim=45)  # Better default viewing angle

    # Increase tick label size
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='z', which='major', labelsize=15)

    # Plot data with consistent colors
    present_patients = np.unique(y)
    for patient in present_patients:
        mask = y == patient
        ax.scatter(X_umap[mask, 0], X_umap[mask, 1], X_umap[mask, 2],
                   color=global_color_map[patient], label=patient, alpha=1,
                   s=50, linewidth=0.3)

    # Axis labels with increased offset
    ax.set_xlabel('UMAP1', labelpad=15, fontweight='bold')
    ax.set_ylabel('UMAP2', labelpad=15, fontweight='bold')
    ax.set_zlabel('UMAP3', labelpad=15, fontweight='bold')
    ax.set_title(f'{cell_type} Networks: 3D UMAP Projection', pad=20, fontsize=14, fontweight='bold')

    # Legend customization
    legend = ax.legend(title="Patient ID", loc='center left',
                       bbox_to_anchor=(1.05, 0.5), frameon=False,
                       title_fontsize=12, fontsize=10)
    legend.get_title().set_fontweight('bold')

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(f"{cell_type}_UMAP_3D.png", dpi=350,
                bbox_inches='tight', facecolor='white')
    plt.close()


# Main execution
if __name__ == '__main__':
    base_dir = os.path.expanduser('~/SingleCellData/LIONESS_Output/')
    cell_types = ['Dendritic', 'Monocyte', 'Progenitor']

    # Create global color scheme first
    all_patients = get_all_patient_ids(base_dir, cell_types)
    color_map = create_global_color_mapping(all_patients)

    # Process each cell type with consistent colors
    for cell_type in cell_types:
        print(f"Generating 3D UMAP for {cell_type}...")
        cell_type_dir = os.path.join(base_dir, cell_type)
        run_umap_3d(cell_type_dir, cell_type, color_map)
