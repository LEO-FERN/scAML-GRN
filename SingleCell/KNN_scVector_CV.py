import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.preprocessing import StandardScaler
from multiprocessing import Pool
from io import StringIO


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


def train_knn(cell_type_dir):
    cell_type = os.path.basename(cell_type_dir)
    print(f"Processing {cell_type} cells...")
    df = process_cell_type(cell_type_dir)

    X = df.drop(['Patient_ID', 'Cell_Type'], axis=1)
    y = df['Patient_ID'].values

    # Initialize cross-validation
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    all_y_test = []
    all_y_pred = []

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        # Split data
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # Train and predict
        knn = KNeighborsClassifier(n_neighbors=3, n_jobs=20)
        knn.fit(X_train_scaled, y_train)
        y_pred = knn.predict(X_test_scaled)

        # Store results
        all_y_test.extend(y_test)
        all_y_pred.extend(y_pred)

        # Fold-specific metrics
        print(f"\nFold {fold_idx + 1} Accuracy: {accuracy_score(y_test, y_pred):.4f}")

    # Generate consolidated reports
    output = StringIO()
    sys.stdout = output

    print("\nConsolidated Classification Report (5-fold CV):")
    print(classification_report(all_y_test, all_y_pred))

    print("\nConfusion Matrix:")
    cm = confusion_matrix(all_y_test, all_y_pred)
    cm_df = pd.DataFrame(cm,
                         index=np.unique(y),
                         columns=np.unique(y))
    cm_df.index.name = 'Actual'
    cm_df.columns.name = 'Predicted'
    print(cm_df)

    sys.stdout = sys.__stdout__

    # Write results
    with open(f"{cell_type}_knn_cv_results.txt", "w") as f:
        f.write(output.getvalue())

    print(f"\nOverall CV Accuracy: {accuracy_score(all_y_test, all_y_pred):.4f}")
    print(f"Results written to {cell_type}_knn_cv_results.txt")


if __name__ == '__main__':
    base_dir = os.path.expanduser('~/SingleCellData/LIONESS_Output/')
    cell_types = ['Dendritic', 'Monocyte', 'Progenitor']

    for cell_type in cell_types:
        cell_type_dir = os.path.join(base_dir, cell_type)
        print(f"\n{'-' * 50}")
        print(f"Processing {cell_type} cells")
        print(f"{'-' * 50}")
        train_knn(cell_type_dir)
