import os
import pandas as pd
import igraph as ig
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Directories and output files
directories = ['~/BinaryFinal/Dendritic', '~/BinaryFinal/Monocyte', '~/BinaryFinal/Progenitor']
#output_files = ['Dendritic_Statistics.csv', 'Monocyte_Statistics.csv', 'Progenitor_Statistics.csv']
output_files = ['Dendritic_Statistics_Consensus.csv', 'Monocyte_Statistics_Consensus.csv', 'Progenitor_Statistics_Consensus.csv']

# Target file keywords
#valid_keywords = ['binary.csv', 'consensus_network.csv']
valid_keywords = ['consensus_network.csv']


def read_network_igraph(csv):
    """Load adjacency matrix and convert to igraph object."""
    df = pd.read_csv(csv, header=0, index_col=False)
    df.index = df.columns  # Ensure index and columns are aligned with the gene names

    # Check if the network is symmetric
    is_symmetric = np.allclose(df.values, df.values.T)

    # Create graph: Nonzero entries represent edges
    g = ig.Graph.Adjacency(
        df.values.tolist(),
        mode=ig.ADJ_UNDIRECTED if is_symmetric else ig.ADJ_DIRECTED
    )
    g.vs["name"] = df.columns  # Assign vertex names (gene names)
    return g, is_symmetric


def compute_statistics_igraph(g, is_symmetric):
    """Compute network statistics using igraph."""
    stats = {}

    # Basic graph metrics
    stats["Nodes"] = g.vcount()
    stats["Edges"] = g.ecount()

    # Components
    if not is_symmetric:
        stats["Weakly Connected Components"] = len(g.clusters(mode="weak"))
        stats["Strongly Connected Components"] = len(g.clusters(mode="strong"))
    else:
        stats["Components"] = len(g.clusters())

    # Largest connected component and diameter
    largest_cc = g.clusters(mode="weak" if not is_symmetric else "strong").giant()
    stats["Diameter"] = largest_cc.diameter()

    # Density and transitivity
    stats["Density"] = g.density()
    stats["Transitivity"] = g.transitivity_undirected() if is_symmetric else g.transitivity_avglocal_undirected()

    # Reciprocity (only for asymmetric graphs)
    stats["Reciprocity"] = g.reciprocity() if not is_symmetric else None

    # Degree statistics
    degrees = g.degree()
    max_degree_idx = np.argmax(degrees)
    max_degree_gene = g.vs[max_degree_idx]["name"] if len(degrees) > 0 else None
    stats.update({
        "Mean Degree": np.mean(degrees),
        "Min Degree": np.min(degrees),
        "Max Degree": np.max(degrees),
        "SD Degree": np.std(degrees),
        "Gene with Highest Degree": max_degree_gene,
    })

    # Clustering coefficient
    clustering_coeffs = g.transitivity_local_undirected(mode="zero")
    stats.update({
        "Mean Clustering Coefficient": np.mean(clustering_coeffs),
        "Variance Clustering Coefficient": np.var(clustering_coeffs),
    })

    # Additional global statistics
    stats["Average Path Length"] = g.average_path_length()
    #stats["Global Efficiency"] = g.global_efficiency()
    stats["Assortativity Coefficient"] = g.assortativity_degree()

    try:
        stats["Girth"] = g.girth()
    except:
        stats["Girth"] = None  # In case the graph is acyclic

    #stats["Clique Number"] = g.clique_number()

    # Modularity (using fast_greedy community detection for undirected graphs)
    #if is_symmetric:
    #    community = g.community_fastgreedy()
    #    stats["Modularity"] = community.modularity
    #else:
    #    stats["Modularity"] = None  # For directed graphs, more complex methods might be needed

    # Edge and Vertex Connectivity
    #stats["Edge Connectivity"] = g.edge_connectivity()
    #stats["Vertex Connectivity"] = g.vertex_connectivity()

    # Number of Triangles
    #stats["Total Triangles"] = sum(g.triangles())

    return stats


def process_file(file_path):
    """Process a single network file and print its stats."""
    g, is_symmetric = read_network_igraph(file_path)
    stats = compute_statistics_igraph(g, is_symmetric)
    stats["Filename"] = os.path.basename(file_path)
    stats["Network Type"] = "Symmetric" if is_symmetric else "Asymmetric"

    # Print stats to follow progress
    print(f"Processed {stats['Filename']} with stats: {stats}")
    return stats


def process_directory_igraph(directory, output_file):
    """Process all valid networks in a directory in parallel and save statistics."""
    expanded_directory = os.path.expanduser(directory)
    valid_files = [
        os.path.join(expanded_directory, file)
        for file in os.listdir(expanded_directory)
        if any(keyword in file for keyword in valid_keywords)
    ]

    network_statistics = []

    with ProcessPoolExecutor(max_workers=30) as executor:
        results = list(executor.map(process_file, valid_files))
        network_statistics.extend(results)

    # Save results to CSV in the same directory as the consensus files
    output_path = os.path.join(expanded_directory, output_file)
    stats_df = pd.DataFrame(network_statistics)
    stats_df.to_csv(output_path, index=False)
    print(f"Statistics saved to {output_path}")


# Main function to sequentially process directories
def main():
    for directory, output_file in zip(directories, output_files):
        process_directory_igraph(directory, output_file)


if __name__ == "__main__":
    main()
