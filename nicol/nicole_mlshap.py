import os
import warnings
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc, classification_report, accuracy_score

from xgboost import XGBClassifier
import shap

warnings.filterwarnings("ignore")


# =========================
# USER SETTINGS
# =========================
INPUT_CSV = r"C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\matched_pcd_probe_labeled_noAmbient.csv"
OUTPUT_DIR = r"C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\xgboost_outputs"

TEST_SIZE = 0.2
RANDOM_STATE = 42

SHAP_EXPLANATION_SIZE = 100
SHAP_TOP_N_FEATURES = 10
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

    numeric_cols = X.select_dtypes(include=[np.number]).columns.tolist()
    X = X[numeric_cols].copy()

    if X.shape[1] == 0:
        raise ValueError("No numeric feature columns found after removing target columns.")

    return X, numeric_cols


def compute_scale_pos_weight(y: pd.Series) -> float:
    """
    For binary classification:
    scale_pos_weight = (# negative) / (# positive)
    """
    y = pd.Series(y)
    neg = (y == 0).sum()
    pos = (y == 1).sum()

    if pos == 0:
        return 1.0
    return float(neg / pos)


def make_pipeline(scale_pos_weight: float) -> Pipeline:
    """
    XGBoost pipeline:
    - median imputation
    - XGBoost classifier

    StandardScaler is not needed for tree-based models.
    """
    model = XGBClassifier(
        n_estimators=300,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=RANDOM_STATE,
        scale_pos_weight=scale_pos_weight if USE_CLASS_WEIGHT_BALANCED else 1.0,
        n_jobs=-1
    )

    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("xgb", model)
    ])
    return pipe


def save_roc_curve(y_true: np.ndarray, y_score: np.ndarray, target_name: str, output_dir: str) -> float:
    """
    Save ROC curve plot and return AUC.
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

    return roc_auc


def save_shap_summary(
    trained_pipeline: Pipeline,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    target_name: str,
    output_dir: str,
) -> pd.DataFrame:
    """
    Save SHAP summary plot for XGBoost.
    Also save top-N mean absolute SHAP feature importances as CSV.
    """
    # Transform data through imputer so SHAP sees the actual model inputs
    X_train_trans = trained_pipeline.named_steps["imputer"].transform(X_train)
    X_test_trans = trained_pipeline.named_steps["imputer"].transform(X_test)

    feature_names = X_train.columns.tolist()

    ex_n = min(SHAP_EXPLANATION_SIZE, X_test_trans.shape[0])
    rng = np.random.default_rng(RANDOM_STATE)
    ex_idx = rng.choice(X_test_trans.shape[0], size=ex_n, replace=False)
    explain_data = X_test_trans[ex_idx]

    xgb_model = trained_pipeline.named_steps["xgb"]

    # TreeExplainer is much faster for XGBoost than KernelExplainer
    explainer = shap.TreeExplainer(xgb_model)
    shap_values = explainer.shap_values(explain_data)

    shap_values = np.asarray(shap_values)

    # For binary XGBoost this is typically (n_samples, n_features)
    if shap_values.ndim != 2:
        raise ValueError(f"Unexpected SHAP output shape: {shap_values.shape}")

    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    shap_importance_df = pd.DataFrame({
        "feature": feature_names,
        "mean_abs_shap": mean_abs_shap
    }).sort_values("mean_abs_shap", ascending=False)

    top_shap_df = shap_importance_df.head(SHAP_TOP_N_FEATURES).reset_index(drop=True)

    csv_path = os.path.join(
        output_dir,
        f"shap_top_{SHAP_TOP_N_FEATURES}_{sanitize_filename(target_name)}.csv"
    )
    top_shap_df.to_csv(csv_path, index=False)

    plt.figure()
    shap.summary_plot(
        shap_values,
        explain_data,
        feature_names=feature_names,
        max_display=SHAP_TOP_N_FEATURES,
        show=False
    )
    plt.title(f"SHAP Summary - Top {SHAP_TOP_N_FEATURES} Features - {target_name}")
    plt.tight_layout()

    out_path = os.path.join(
        output_dir,
        f"shap_summary_top_{SHAP_TOP_N_FEATURES}_{sanitize_filename(target_name)}.png"
    )
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    return top_shap_df


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

        unique_classes = sorted(pd.Series(y_valid).dropna().unique().tolist())
        if len(unique_classes) != 2:
            print(f"Skipping {target_col}: target is not binary. Found classes: {unique_classes}")
            continue

        class_counts = y_valid.value_counts().sort_index()
        print(f"Class counts for {target_col}:")
        print(class_counts)

        X_train, X_test, y_train, y_test = train_test_split(
            X_valid,
            y_valid,
            test_size=TEST_SIZE,
            random_state=RANDOM_STATE,
            stratify=y_valid
        )

        scale_pos_weight = compute_scale_pos_weight(y_train)
        print(f"scale_pos_weight for {target_col}: {scale_pos_weight:.4f}")

        pipeline = make_pipeline(scale_pos_weight=scale_pos_weight)
        pipeline.fit(X_train, y_train)

        y_proba = pipeline.predict_proba(X_test)[:, 1]
        y_pred = pipeline.predict(X_test)

        roc_auc = save_roc_curve(y_test.to_numpy(), y_proba, target_col, OUTPUT_DIR)
        acc = accuracy_score(y_test, y_pred)

        print(f"AUC for {target_col}: {roc_auc:.4f}")
        print(f"Accuracy for {target_col}: {acc:.4f}")
        print(classification_report(y_test, y_pred, digits=4))

        try:
            top_shap_df = save_shap_summary(pipeline, X_train, X_test, target_col, OUTPUT_DIR)
            shap_status = "saved"

            print(f"\nTop {SHAP_TOP_N_FEATURES} SHAP features for {target_col}:")
            print(top_shap_df.to_string(index=False))

            top_features_str = "; ".join(top_shap_df["feature"].tolist())
        except Exception as e:
            print(f"SHAP failed for {target_col}: {e}")
            shap_status = f"failed: {e}"
            top_features_str = ""

        summary_rows.append({
            "target": target_col,
            "n_features": X_train.shape[1],
            "n_train": X_train.shape[0],
            "n_test": X_test.shape[0],
            "auc": roc_auc,
            "accuracy": acc,
            "scale_pos_weight": scale_pos_weight,
            "shap_status": shap_status,
            "top_features": top_features_str
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTPUT_DIR, "model_summary.csv")
    summary_df.to_csv(summary_path, index=False)

    print("\nDone.")
    print(f"Outputs saved in: {OUTPUT_DIR}")
    print(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()