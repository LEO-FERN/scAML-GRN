import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

def parse_confusion_matrix(file_path):
    """Parse confusion matrix from results file with robust header handling"""
    with open(file_path, 'r') as f:
        lines = [line.rstrip() for line in f]

    try:
        cm_start = next(i for i, line in enumerate(lines) if line.startswith("Confusion Matrix:"))
    except StopIteration:
        raise ValueError(f"Confusion Matrix section missing in {file_path}")

    header_line = None
    for i in range(cm_start + 1, cm_start + 4):
        if i < len(lines) and "Predicted" in lines[i]:
            header_line = lines[i]
            break

    if not header_line:
        raise ValueError(f"Predicted header not found in {file_path}")

    predicted_classes = re.split(r'\s{2,}', header_line.split("Predicted")[-1].strip())

    matrix_data = []
    actual_classes = []
    data_start = cm_start + 3

    for line in lines[data_start:]:
        if not line.strip() or line.startswith("Average"):
            break
        parts = re.split(r'\s{2,}', line.strip())
        actual_classes.append(parts[0])
        matrix_data.append([int(x) for x in parts[1:len(predicted_classes) + 1]])

    return pd.DataFrame(matrix_data, index=actual_classes, columns=predicted_classes)

def create_heatmap(ax, df, title, show_cbar=False, cbar_ax=None):
    """Create a single heatmap within a subplot"""
    df_percentages = df.div(df.sum(axis=1), axis=0) * 100

    colors = ['#FFFFFF', '#FFF0F0', '#FFE0E0', '#FFD0D0', '#FFC0C0', '#FFB0B0',
              '#FFA0A0', '#FF9090', '#FF8080', '#FF7070', '#FF6060', '#FF5050',
              '#FF4040', '#FF3030', '#FF2020', '#FF1010', '#FF0000']
    n_bins = len(colors)
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_red", colors, N=n_bins)

    sns.heatmap(
        df_percentages,
        annot=False,
        fmt=".1f",
        cmap=cmap,
        cbar=show_cbar,
        cbar_ax=cbar_ax,
        linewidths=0.5,
        linecolor='black',
        vmin=0,
        vmax=100,
        ax=ax,
        square=True  # Ensures square heatmap
    )

    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel("")  # Remove x label
    ax.set_ylabel("")  # Remove y label
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)  # Reduced font size
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)  # Reduced font size
    ax.yaxis.set_tick_params(pad=1)  # Increased padding

def create_combined_heatmap(cell_type, heatmap_data, output_path):
    """Create a combined heatmap image with KNN, SVM, and RF side by side, one legend at the end (equal sizes)."""
    # 3 heatmaps + 1 colorbar column
    fig = plt.figure(figsize=(20, 6))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.06], wspace=0.25)

    axes = [plt.subplot(gs[i]) for i in range(3)]
    cbar_ax = plt.subplot(gs[3])

    ordered_classifiers = ["KNN", "SVM", "RF"]
    for idx, (ax, classifier) in enumerate(zip(axes, ordered_classifiers)):
        # Only show colorbar on last plot, and direct it to the dedicated cbar_ax
        show_cbar = (idx == 2)
        create_heatmap(
            ax,
            heatmap_data[classifier],
            f"{cell_type} - {classifier}",
            show_cbar=True if show_cbar else False,
            cbar_ax=cbar_ax if show_cbar else None
        )

    cbar_ax.set_ylabel("Percentage (%)", fontsize=12)
    cbar_ax.yaxis.set_label_position('right')
    cbar_ax.yaxis.tick_right()

    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='svg')
    plt.close()

def main():
    output_dir = Path("heatmaps")
    output_dir.mkdir(exist_ok=True)

    cell_types = {"Dendritic": {}, "Monocyte": {}, "Progenitor": {}}
    classifiers = {"knn": "KNN", "svm": "SVM", "rf": "RF"}

    for file in os.listdir():
        if file.endswith("_cv_results.txt"):
            parts = file.split("_")
            cell_type = parts[0]
            classifier = parts[1].lower()

            if cell_type in cell_types and classifier in classifiers:
                cell_types[cell_type][classifiers[classifier]] = parse_confusion_matrix(file)

    for cell_type, heatmap_data in cell_types.items():
        if set(heatmap_data.keys()) == {"KNN", "SVM", "RF"}:  # Ensure all three classifiers are available
            output_path = output_dir / f"{cell_type}_combined_heatmap.svg"
            create_combined_heatmap(cell_type, heatmap_data, output_path)
            print(f"Created: {output_path}")

if __name__ == "__main__":
    main()
