import pandas as pd
import numpy as np
import os
import re
import argparse

def get_network_IDs(path):
    files = os.listdir(path)

    IDs = []
    for file in files:
        AML_match = re.search(r"AML[^_]+_", file)
        BM_match = re.search(r"BM[^_]+_", file)
        if AML_match:
            IDs.append(AML_match.group(0).strip("_"))
        if BM_match:
            IDs.append(BM_match.group(0).strip("_"))

    return files, list(set(IDs))

def network_focus(path):
    files, IDs = get_network_IDs(path)
    focus = []

    for ID in IDs:
        matching_files = [file for file in files if ID in file]
        # Drop files that contain 'matrix.csv'
        matching_files = [file for file in matching_files if 'matrix.csv' not in file]
        matching_files = [file for file in matching_files if '.csv' in file]
        matching_files.append(ID)

        focus.append(matching_files)
    return focus

def read_network(file_path):
    """Read network data from a CSV file."""
    df = pd.read_csv(file_path, header=0, index_col=False)
    df.index = df.columns  # Ensure index and columns are aligned with the gene names
    return df

def apply_threshold(network, threshold):
    """Apply threshold to create a binary network."""
    return (network > threshold).astype(int)

def save_binary_network(network, output_path):
    """Save the binary network to a CSV file."""
    network.to_csv(output_path, index=False)

def process_networks(input_dir, output_dir, thresholds):
    """Process all networks, apply thresholds, and save binary networks."""
    #methods = ['ARACNE.', 'CLR.', 'MRNET.', 'GENIE.']
    #methods = ['ARACNE.', 'CLR.', 'MRNET.', 'GENIE_SYM.']
    methods = ['GENIE_SYM.']

    focus_list = network_focus(input_dir)

    for focus in focus_list:
        ID = focus[-1]  # The last element is the ID
        for file in focus[:-1]:  # Exclude the last element (ID) when processing files
            for method in methods:
                if method in file:
                    input_file = os.path.join(input_dir, file)
                    output_file = os.path.join(output_dir, f"{ID}_{method.strip('.')}_binary.csv")

                    # Read the network
                    network = read_network(input_file)

                    # Apply threshold
                    binary_network = apply_threshold(network, thresholds[method.strip('.')])

                    # Save binary network
                    save_binary_network(binary_network, output_file)

                    print(f"Processed {ID} {method.strip('.')} network")

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate binary networks for specific IDs.")
    parser.add_argument("--input_dir", type=str, help="Input directory containing network files")
    parser.add_argument("--output_dir", type=str, help="Output directory for binary networks")
    args = parser.parse_args()

    input_directory = os.path.expanduser(args.input_dir)
    output_directory = os.path.expanduser(args.output_dir)

    # Set your thresholds here
    thresholds = {
        'ARACNE': 0,
        'CLR': 2.84166785553267,
        'MRNET': 0.00306545740959235,
        'GENIE_SYM': 0.000288810502691997
    }

    # Ensure output directory exists
    #os.makedirs(output_directory, exist_ok=True)

    # Process networks
    process_networks(input_directory, output_directory, thresholds)

    print("All networks processed and binary networks saved.")