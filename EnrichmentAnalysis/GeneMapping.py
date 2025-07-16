import numpy as np
import pandas as pd
import csv
import os

def create_gene_index_mapping(gene_list_file, output_csv):
    # Read the gene list
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f]

    # Calculate the size of the network
    n = len(genes)

    # Create the mapping
    mapping = []
    for i in range(n):
        for j in range(i+1, n):
            index = len(mapping)
            mapping.append((index, genes[i], genes[j]))

    # Write to CSV
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Index', 'Gene 1', 'Gene 2'])
        writer.writerows(mapping)

    print(f"Mapping saved to {output_csv}")
    print(f"Total number of features: {len(mapping)}")

# Usage - replace monocyte, progenitor, dendritic etc...
gene_list_file = os.path.expanduser('~/EnrichmentData/progenitor_genes.txt')
output_csv = os.path.expanduser('~/EnrichmentData/progenitor_gene_index_mapping.csv')

create_gene_index_mapping(gene_list_file, output_csv)
