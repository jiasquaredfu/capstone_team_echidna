import os
import warnings
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc, classification_report

import shap

warnings.filterwarnings("ignore")


# =========================
# USER SETTINGS
# =========================
INPUT_CSV = r"C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\matched_pcd_probe_labeled_noAmbient.csv"
OUTPUT_DIR = r"C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\svm_outputs"

TEST_SIZE = 0.2
RANDOM_STATE = 42

# SHAP can be slow for nonlinear SVMs, so keep these modest
SHAP_BACKGROUND_SIZE = 50
SHAP_EXPLANATION_SIZE = 50

# If your classes are imbalanced, this usually helps
USE_CLASS_WEIGHT_BALANCED = True


# =========================
# HELPERS
# =========================
def sanitize_filename(name: str) -> str:
    """Make a safe filename from a column name."""
    bad_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*', ' ']
    out = name
    for ch in bad_chars:
        out = out.replace(ch, "_")
    return out


def get_target_columns(df: pd.DataFrame) -> List[str]:
    """
    Use the last 3 columns as targets.
    """
    if df.shape[1] < 4:
        raise ValueError("The CSV does not have enough columns to split features and 3 targets.")
    return list(df.columns[-3:])


def build_feature_matrix(
    df: pd.DataFrame,
    target_cols: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build feature set by removing the 3 target columns and keeping only numeric features.
    """
    X = df.drop(columns=target_cols).copy()

    # Keep numeric columns only
    numeric_cols = X.select_dtypes(include=[np.number]).columns.tolist()
    X = X[numeric_cols].copy()

    if X.shape[1] == 0:
        raise ValueError("No numeric feature columns found after removing target columns.")

    return X, numeric_cols


def make_pipeline() -> Pipeline:
    """
    Standard SVM pipeline:
    - median imputation
    - standardization
    - RBF SVM with probabilities enabled for ROC + SHAP
    """
    class_weight = "balanced" if USE_CLASS_WEIGHT_BALANCED else None

    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("svm", SVC(
            kernel="rbf",
            probability=True,
            class_weight=class_weight,
            random_state=RANDOM_STATE
        ))
    ])
    return pipe


def save_roc_curve(y_true: np.ndarray, y_score: np.ndarray, target_name: str, output_dir: str) -> None:
    """
    Save ROC curve plot.
    """
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.3f}")
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve - {target_name}")
    plt.legend(loc="lower right")
    plt.tight_layout()

    out_path = os.path.join(output_dir, f"roc_{sanitize_filename(target_name)}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def save_shap_summary(
    trained_pipeline: Pipeline,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    target_name: str,
    output_dir: str,
) -> None:
    """
    Save SHAP summary plot for nonlinear SVM using KernelExplainer.

    Since KernelExplainer can be slow, this uses small samples from train/test.
    """
    # Transform data through imputer + scaler so SHAP explains the model input space
    X_train_trans = trained_pipeline[:-1].transform(X_train)
    X_test_trans = trained_pipeline[:-1].transform(X_test)

    feature_names = X_train.columns.tolist()

    # Sample background and explanation sets
    rng = np.random.default_rng(RANDOM_STATE)

    bg_n = min(SHAP_BACKGROUND_SIZE, X_train_trans.shape[0])
    ex_n = min(SHAP_EXPLANATION_SIZE, X_test_trans.shape[0])

    bg_idx = rng.choice(X_train_trans.shape[0], size=bg_n, replace=False)
    ex_idx = rng.choice(X_test_trans.shape[0], size=ex_n, replace=False)

    background = X_train_trans[bg_idx]
    explain_data = X_test_trans[ex_idx]

    svm_model = trained_pipeline.named_steps["svm"]

    # Explain probability of class 1
    explainer = shap.KernelExplainer(
        svm_model.predict_proba,
        background
    )

    shap_values = explainer.shap_values(explain_data)

    # For binary classification, SHAP may return:
    # - a list [class0, class1]
    # or
    # - a single array depending on version
    if isinstance(shap_values, list):
        shap_for_positive = shap_values[1]
    else:
        # Some SHAP versions return array directly
        shap_for_positive = shap_values

    plt.figure()
    shap.summary_plot(
        shap_for_positive,
        explain_data,
        feature_names=feature_names,
        show=False
    )
    plt.title(f"SHAP Summary - {target_name}")
    plt.tight_layout()

    out_path = os.path.join(output_dir, f"shap_summary_{sanitize_filename(target_name)}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


# =========================
# MAIN
# =========================
def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)

    # Last 3 columns are targets
    target_cols = get_target_columns(df)
    print(f"Target columns: {target_cols}")

    # Build features once, excluding all target columns
    X, feature_names = build_feature_matrix(df, target_cols)
    print(f"Number of feature columns used: {len(feature_names)}")

    summary_rows = []

    for target_col in target_cols:
        print("\n" + "=" * 70)
        print(f"Training model for target: {target_col}")

        y = df[target_col].copy()

        # Drop rows with missing target
        valid_mask = ~pd.isna(y)
        X_valid = X.loc[valid_mask].copy()
        y_valid = y.loc[valid_mask].astype(int).copy()

        # Basic binary check
        unique_classes = sorted(pd.Series(y_valid).dropna().unique().tolist())
        if len(unique_classes) != 2:
            print(f"Skipping {target_col}: target is not binary. Found classes: {unique_classes}")
            continue

        # Stratified split
        X_train, X_test, y_train, y_test = train_test_split(
            X_valid,
            y_valid,
            test_size=TEST_SIZE,
            random_state=RANDOM_STATE,
            stratify=y_valid
        )

        # Train
        pipeline = make_pipeline()
        pipeline.fit(X_train, y_train)

        # Predict probabilities for ROC
        y_proba = pipeline.predict_proba(X_test)[:, 1]
        y_pred = pipeline.predict(X_test)

        # ROC
        fpr, tpr, _ = roc_curve(y_test, y_proba)
        roc_auc = auc(fpr, tpr)

        print(f"AUC for {target_col}: {roc_auc:.4f}")
        print(classification_report(y_test, y_pred, digits=4))

        # Save ROC
        save_roc_curve(y_test.to_numpy(), y_proba, target_col, OUTPUT_DIR)

        # Save SHAP
        try:
            save_shap_summary(pipeline, X_train, X_test, target_col, OUTPUT_DIR)
            shap_status = "saved"
        except Exception as e:
            print(f"SHAP failed for {target_col}: {e}")
            shap_status = f"failed: {e}"

        summary_rows.append({
            "target": target_col,
            "n_features": X_train.shape[1],
            "n_train": X_train.shape[0],
            "n_test": X_test.shape[0],
            "auc": roc_auc,
            "shap_status": shap_status
        })

    # Save summary
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTPUT_DIR, "model_summary.csv")
    summary_df.to_csv(summary_path, index=False)

    print("\nDone.")
    print(f"Outputs saved in: {OUTPUT_DIR}")
    print(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()