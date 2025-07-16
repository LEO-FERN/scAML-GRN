import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
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

def train_svm(cell_type_dir):
    cell_type = os.path.basename(cell_type_dir)
    print(f"Processing {cell_type} cells...")
    df = process_cell_type(cell_type_dir)

    X = df.drop(['Patient_ID', 'Cell_Type'], axis=1).values
    y = df['Patient_ID'].values

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    all_y_test = []
    all_y_pred = []
    accuracies = []

    for fold_idx, (train_index, test_index) in enumerate(skf.split(X, y), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        svm_model = SVC(kernel='rbf', random_state=42)
        svm_model.fit(X_train_scaled, y_train)
        y_pred = svm_model.predict(X_test_scaled)

        all_y_test.extend(y_test)
        all_y_pred.extend(y_pred)
        fold_accuracy = svm_model.score(X_test_scaled, y_test)
        accuracies.append(fold_accuracy)
        print(f"Fold {fold_idx} Accuracy: {fold_accuracy:.4f}")

    output = StringIO()
    sys.stdout = output

    print("\nClassification Report (5-fold CV):")
    print(classification_report(all_y_test, all_y_pred))

    cm = confusion_matrix(all_y_test, all_y_pred)
    patient_ids = np.unique(y)
    cm_df = pd.DataFrame(cm, index=patient_ids, columns=patient_ids)
    cm_df.index.name = 'Actual'
    cm_df.columns.name = 'Predicted'

    print("\nConfusion Matrix:")
    print(cm_df)
    print(f"\nAverage Accuracy: {np.mean(accuracies):.4f} ± {np.std(accuracies):.4f}")

    sys.stdout = sys.__stdout__

    with open(f"{cell_type}_svm_cv_results.txt", "w") as f:
        f.write(output.getvalue())

    print(f"Results written to {cell_type}_svm_cv_results.txt")

if __name__ == '__main__':
    base_dir = os.path.expanduser('~/SingleCellData/LIONESS_Output/')
    cell_types = ['Dendritic', 'Monocyte', 'Progenitor']

    for cell_type in cell_types:
        cell_type_dir = os.path.join(base_dir, cell_type)
        print(f"\nProcessing {cell_type} cells")
        print("=" * 50)
        train_svm(cell_type_dir)
        print("\n")
