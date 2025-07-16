import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Read the CSV files
dendritic = pd.read_csv('~/BinaryFinal/Dendritic/Dendritic_Statistics_Consensus.csv')
progenitor = pd.read_csv('~/BinaryFinal/Progenitor/Progenitor_Statistics_Consensus.csv')
monocyte = pd.read_csv('~/BinaryFinal/Monocyte/Monocyte_Statistics_Consensus.csv')

# Combine the dataframes and add a 'Cell_Type' column
dendritic['Cell_Type'] = 'Dendritic'
progenitor['Cell_Type'] = 'Progenitor'
monocyte['Cell_Type'] = 'Monocyte'
combined_df = pd.concat([dendritic, progenitor, monocyte], ignore_index=True)

# Extract network IDs from the Filename column
combined_df['Network_ID'] = combined_df['Filename'].apply(lambda x: re.search(r'(AML\d+[A-Z]?|BM\d+)', x).group(1))

# Select numerical features for PCA
numerical_features = ['Edges', 'Density', 'Transitivity', 'Mean Degree', 'Min Degree', 'Max Degree',
                      'SD Degree', 'Mean Clustering Coefficient', 'Variance Clustering Coefficient',
                      'Average Path Length', 'Assortativity Coefficient'] # 'Nodes', 'Diameter', 'Reciprocity', 'Girth'

# Prepare the data for PCA
X = combined_df[numerical_features]
y = combined_df['Cell_Type']

# Standardize the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Create a dataframe with PCA results
pca_df = pd.DataFrame(data=X_pca, columns=['PC1', 'PC2'])
pca_df['Cell_Type'] = y
pca_df['Network_ID'] = combined_df['Network_ID']

# Create plot with Patient ID color mapping
plt.figure(figsize=(12, 8))

# Get unique patient IDs and assign colors
#patient_ids = sorted(pca_df['Network_ID'].unique())
#colors = plt.cm.tab20(np.linspace(0, 1, len(patient_ids)))
#color_map = dict(zip(patient_ids, colors))

# Create color mapping for Network_IDs
patient_ids = sorted(pca_df['Network_ID'].unique())
colors = plt.cm.get_cmap('tab20', len(patient_ids))
color_map = {pid: colors(i) for i, pid in enumerate(patient_ids)}

# Define markers for cell types (but won't appear in legend)
markers = {'Dendritic': 'o', 'Progenitor': 's', 'Monocyte': '^'}

# Plot each cell type with its marker, colored by patient ID
for cell_type, marker in markers.items():
    subset = pca_df[pca_df['Cell_Type'] == cell_type]
    plt.scatter(subset['PC1'], subset['PC2'],
                c=subset['Network_ID'].map(color_map),
                marker=marker,
                s=150,
                edgecolor='w',
                linewidth=0.5,
                label='_nolegend_')  # This hides these from legend

# Create legend only for Patient IDs
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=pid,
                             markerfacecolor=color, markersize=15)
                  for pid, color in color_map.items()]

plt.legend(handles=legend_elements, title='Patient IDs',
           loc='center left', bbox_to_anchor=(1, 0.5))

plt.title('PCA of Network Statistics')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pca_statistics_enhanced.png', dpi=300, bbox_inches='tight')