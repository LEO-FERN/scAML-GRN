import os
import pandas as pd
import csv


def expand_path(path):
    """Expand ~ to full home directory path"""
    return os.path.expanduser(path)


def process_networks(cell_dirs, output_dir):
    """
    Process multiple cell type directories and save results to specified output

    Args:
        cell_dirs (list): List of paths to cell type directories
        output_dir (str): Destination directory for output file
    """
    # Expand and create output directory if needed
    output_dir = expand_path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'top_connected_genes.csv')

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['CellType', 'PatientID', 'TopGenes'])

        for cell_dir in cell_dirs:
            # Expand cell directory path
            cell_dir = expand_path(cell_dir)
            cell_type = os.path.basename(os.path.normpath(cell_dir))

            for file in os.listdir(cell_dir):
                if file.endswith('_consensus_network.csv'):
                    patient_id = file.split('_consensus_network')[0]

                    file_path = os.path.join(cell_dir, file)
                    adj_matrix = pd.read_csv(file_path, header=0, index_col=None)
                    adj_matrix.index = adj_matrix.columns

                    connectivity = adj_matrix.sum(axis=1)
                    top_genes = connectivity.nlargest(2500).index.tolist()

                    writer.writerow([
                        cell_type,
                        patient_id,
                        ','.join(top_genes)
                    ])


if __name__ == '__main__':

    print('Top Genes...')

    CELL_TYPE_DIRS = [
        '~/BinaryFinal/Dendritic',
        '~/BinaryFinal/Monocyte',
        '~/BinaryFinal/Progenitor'
    ]

    OUTPUT_DIR = '~/EnrichmentData'

    process_networks(
        cell_dirs=CELL_TYPE_DIRS,
        output_dir=OUTPUT_DIR
    )
