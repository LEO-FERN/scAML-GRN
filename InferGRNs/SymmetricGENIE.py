import os
import numpy as np
import pandas as pd
import glob
from multiprocessing import Pool, cpu_count

def process_file(file):
    # Read the CSV file
    network = pd.read_csv(file, header=0, index_col=False)
    network.index = network.columns

    # Convert to numpy array
    network_array = network.to_numpy()

    # Make symmetric
    symmetric_network = (network_array + network_array.T) / 2

    # Convert back to DataFrame
    symmetric_df = pd.DataFrame(symmetric_network, index=network.index, columns=network.columns)

    # Create output filename
    output_file = file.replace("_GENIE.csv", "_GENIE_SYM.csv")

    # Save to CSV
    symmetric_df.to_csv(output_file, index=False)

    return f"Processed {file} -> {output_file}"

def main():
    # Define the directories
    directories = [
        os.path.expanduser("~/Data/Final_Dendritic_Net"),
        os.path.expanduser("~/Data/Final_Monocyte_Net"),
        os.path.expanduser("~/Data/Final_Progenitor_Net")
    ]

    # Collect all input files
    all_input_files = []
    for directory in directories:
        input_files = glob.glob(os.path.join(directory, "*_GENIE.csv"))
        all_input_files.extend(input_files)

    # Determine the number of cores to use (max 8)
    num_cores = min(8, cpu_count())

    # Create a pool of workers
    with Pool(num_cores) as pool:
        # Map the process_file function to all input files
        results = pool.map(process_file, all_input_files)

    # Print results
    for result in results:
        print(result)

    print("All files processed successfully.")

if __name__ == "__main__":
    main()
