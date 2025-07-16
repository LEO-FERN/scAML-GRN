import numpy as np
import pandas as pd
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
        matching_files = [file for file in matching_files if 'binary.csv' in file]
        matching_files.append(ID)

        focus.append(matching_files)
    return focus

def get_consensus(path):
    focus = network_focus(path)

    for networks in focus:
        # Load the networks
        net1 = pd.read_csv(f'{path}/{networks[0]}', header=0, index_col=None)
        net1.index = net1.columns
        net2 = pd.read_csv(f'{path}/{networks[1]}', header=0, index_col=None)
        net2.index = net2.columns
        net3 = pd.read_csv(f'{path}/{networks[2]}', header=0, index_col=None)
        net3.index = net3.columns
        net4 = pd.read_csv(f'{path}/{networks[3]}', header=0, index_col=None)
        net4.index = net4.columns
        ID = networks[-1]

        # Extract the gene names
        genes = net1.index

        # Convert to numpy arrays
        matrix1 = net1.to_numpy(dtype=bool)
        matrix2 = net2.to_numpy(dtype=bool)
        matrix3 = net3.to_numpy(dtype=bool)
        matrix4 = net4.to_numpy(dtype=bool)

        # Create the consensus network (union of all edges)
        consensus_matrix = np.logical_or.reduce([matrix1, matrix2, matrix3, matrix4])

        # Convert the consensus network back to a DataFrame
        consensus_df = pd.DataFrame(consensus_matrix.astype(int), columns=genes, index=genes)

        # Save the consensus network
        consensus_df.to_csv(f"{path}/{ID}_consensus_network.csv", index = False)
        print(f'{ID}_consensus_network.csv created')

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate consensus networks from input files.")
    parser.add_argument("path", type=str, help="Path to the directory containing network files.")
    args = parser.parse_args()

    # Validate path
    if not os.path.isdir(args.path):
        print(f"Error: The provided path '{args.path}' is not a valid directory.")
        return

    # Run the consensus network generation
    print(f"Processing networks in directory: {args.path}")
    get_consensus(args.path)

if __name__ == "__main__":
    main()
