#!/usr/bin/env python3
"""
train_crispr_model.py

Train a gradient boosting CRISPR read classifier from the labeled feature dataset
produced by build_crispr_training_dataset.py.

Cross-validation strategy: leave-one-individual-out (LOIO) or leave-one-dataset-out
to prevent data leakage between biological replicates.

Output: assets/models/crispr_model_v*.pkl  (NEVER overwrites v1)

Usage:
  python train_crispr_model.py \\
    --training-data training_data/crispr_training_v2.parquet \\
    --output-model assets/models/crispr_model_v2.pkl \\
    --cv-by individual_id \\
    --threshold 0.70
"""

import argparse
import os
import sys
import warnings

# Add repo root to sys.path so locally-installed packages (pyarrow, joblib, xgboost, etc.) are found
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import (
    roc_auc_score, precision_score, recall_score, f1_score,
    confusion_matrix, classification_report
)
from sklearn.preprocessing import label_binarize

warnings.filterwarnings('ignore', category=UserWarning)

FEATURE_NAMES = [
    'read_pair_gap', 'read_softclip',
    'read_del_vs_control', 'read_ins_vs_control', 'read_mismatch_vs_control',
    'deletion_exclusive_to_edited',
    'is_at_any_target_site', 'distance_to_closest_pam',
    'total_indel_size', 'indel_size_category', 'indel_complexity_score',
    'softclip_dist_to_PAM', 'is_in_homopolymer', 'dist_to_microsatellite',
]

# ---------------------------------------------------------------------------
# Model factory — prefer XGBoost, fall back to sklearn GradientBoosting
# ---------------------------------------------------------------------------

def _make_model(n_estimators=500, max_depth=6, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8, random_state=42,
                scale_pos_weight=1.0):
    try:
        from xgboost import XGBClassifier
        return XGBClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            subsample=subsample,
            colsample_bytree=colsample_bytree,
            scale_pos_weight=scale_pos_weight,
            eval_metric='logloss',
            use_label_encoder=False,
            random_state=random_state,
            n_jobs=-1,
            verbosity=0,
        )
    except ImportError:
        from sklearn.ensemble import GradientBoostingClassifier
        print("XGBoost not available — using sklearn GradientBoostingClassifier", flush=True)
        return GradientBoostingClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            subsample=subsample,
            random_state=random_state,
            verbose=0,
        )


# ---------------------------------------------------------------------------
# Leave-one-group-out cross-validation
# ---------------------------------------------------------------------------

def run_cv(df: pd.DataFrame, cv_col: str, threshold: float, neg_pos_ratio=None) -> dict:
    groups = df[cv_col].unique()
    print(f"\nLeave-one-{cv_col}-out CV: {len(groups)} folds", flush=True)

    fold_results = []
    all_true = []
    all_proba = []

    for g in sorted(groups):
        train_mask = df[cv_col] != g
        test_mask  = df[cv_col] == g

        X_train = df.loc[train_mask, FEATURE_NAMES].values
        y_train = df.loc[train_mask, 'label'].values
        X_test  = df.loc[test_mask,  FEATURE_NAMES].values
        y_test  = df.loc[test_mask,  'label'].values

        if len(np.unique(y_test)) < 2:
            print(f"  Fold {g}: only one class in test set — skipping", flush=True)
            continue

        n_pos = y_train.sum()
        n_neg = len(y_train) - n_pos
        # If negatives were already undersampled, scale_pos_weight is close to
        # neg_pos_ratio and should not be doubled up; set it explicitly.
        spw = (n_neg / max(n_pos, 1)) if neg_pos_ratio is None else 1.0

        model = _make_model(scale_pos_weight=spw)
        model.fit(X_train, y_train)

        proba  = model.predict_proba(X_test)[:, 1]
        pred   = (proba >= threshold).astype(int)

        auc    = roc_auc_score(y_test, proba)
        prec   = precision_score(y_test, pred, zero_division=0)
        rec    = recall_score(y_test, pred, zero_division=0)
        f1     = f1_score(y_test, pred, zero_division=0)
        cm     = confusion_matrix(y_test, pred)

        fold_results.append({'fold': g, 'auc': auc, 'precision': prec,
                              'recall': rec, 'f1': f1, 'n_test': len(y_test),
                              'n_pos_test': y_test.sum()})
        all_true.extend(y_test)
        all_proba.extend(proba)

        print(f"  {g:<30s}  AUC={auc:.3f}  P={prec:.3f}  R={rec:.3f}  F1={f1:.3f}  "
              f"n={len(y_test):,} ({y_test.sum()} pos)  CM={cm.tolist()}", flush=True)

    fold_df = pd.DataFrame(fold_results)
    print(f"\nCV Summary (mean ± std):", flush=True)
    for col in ('auc', 'precision', 'recall', 'f1'):
        vals = fold_df[col]
        print(f"  {col:<12s}: {vals.mean():.3f} ± {vals.std():.3f}", flush=True)

    overall_auc = roc_auc_score(all_true, all_proba)
    print(f"\n  Overall AUC (pooled): {overall_auc:.4f}", flush=True)

    return {
        'fold_results': fold_df,
        'overall_auc': overall_auc,
        'all_true': np.array(all_true),
        'all_proba': np.array(all_proba),
    }


# ---------------------------------------------------------------------------
# Feature importance
# ---------------------------------------------------------------------------

def print_feature_importance(model, top_n=17):
    if hasattr(model, 'feature_importances_'):
        imp = pd.Series(model.feature_importances_, index=FEATURE_NAMES)
        imp = imp.sort_values(ascending=False)
        print("\nFeature importances (top {top_n}):".format(top_n=top_n))
        for feat, val in imp.head(top_n).items():
            bar = '#' * int(val * 60)
            print(f"  {feat:<40s} {val:.4f}  {bar}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Train CRISPR read classifier (v2)")
    p.add_argument('--training-data',  required=True, nargs='+',
                   help='Input parquet file(s) from build_crispr_training_dataset.py '
                        '(accepts multiple files, which are concatenated)')
    p.add_argument('--output-model',   required=True,
                   help='Output model path (e.g. assets/models/crispr_model_v2.pkl)')
    p.add_argument('--cv-by',          default='individual_id',
                   choices=['individual_id', 'dataset'],
                   help='Cross-validation grouping column (default: individual_id)')
    p.add_argument('--threshold',      type=float, default=0.70,
                   help='Decision threshold for binary classification (default: 0.70)')
    p.add_argument('--n-estimators',   type=int,   default=500)
    p.add_argument('--max-depth',      type=int,   default=6)
    p.add_argument('--learning-rate',  type=float, default=0.05)
    p.add_argument('--neg-pos-ratio',  type=float, default=None,
                   help='Undersample negatives to this multiple of positives, '
                        'stratified by --cv-by column. E.g. 10 → 10 negatives per '
                        'positive. Default: use all negatives (with scale_pos_weight).')
    p.add_argument('--no-cv',          action='store_true',
                   help='Skip cross-validation, train final model directly')
    return p.parse_args()


def load_data_balanced(paths, neg_pos_ratio=None, random_state=42):
    """
    Load parquet files without ever holding all rows in memory at once.

    When neg_pos_ratio is set:
      Pass 1 — scan only the 'label' column (fast, tiny RAM) to count
               positives and negatives, compute the uniform sampling fraction.
      Pass 2 — iterate row-group batches; keep every positive row, randomly
               sample negatives at neg_frac per batch.

    Peak RAM ≈ one batch (~200 MB) + accumulated positives (~50 MB for 205k rows)
              + accumulated sampled negatives (~200 MB for 2 M rows at 1:10).
    Total well under 2 GB regardless of source file size.

    When neg_pos_ratio is None all rows are loaded (original behaviour).
    """
    import pyarrow.parquet as pq

    if neg_pos_ratio is not None:
        # Pass 1: count labels only
        n_pos_total = n_neg_total = 0
        for path in paths:
            pf = pq.ParquetFile(path)
            for batch in pf.iter_batches(batch_size=500_000, columns=['label']):
                arr = batch.column('label').to_pylist()
                n_pos_total += arr.count(1)
                n_neg_total += arr.count(0)
        n_target  = int(n_pos_total * neg_pos_ratio)
        neg_frac  = min(1.0, n_target / max(n_neg_total, 1))
        print(f"  Source: {n_pos_total:,} positives, {n_neg_total:,} negatives", flush=True)
        print(f"  Sampling {neg_frac:.4f} of negatives → target ratio 1:{neg_pos_ratio}", flush=True)
    else:
        neg_frac = 1.0

    # Pass 2: load in batches, subsample negatives
    rng = np.random.default_rng(random_state)
    pos_chunks = []
    neg_chunks = []

    for path in paths:
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(batch_size=500_000):
            df_batch = batch.to_pandas()
            pos_chunks.append(df_batch[df_batch['label'] == 1].copy())
            neg_batch = df_batch[df_batch['label'] == 0]
            if neg_frac < 1.0 and len(neg_batch):
                neg_batch = neg_batch.sample(frac=neg_frac, random_state=random_state)
            neg_chunks.append(neg_batch)
            del df_batch

    pos = pd.concat(pos_chunks, ignore_index=True)
    neg = pd.concat(neg_chunks, ignore_index=True)
    del pos_chunks, neg_chunks

    df = pd.concat([pos, neg], ignore_index=True)
    del pos, neg
    df = df.sample(frac=1, random_state=random_state).reset_index(drop=True)
    return df


def main():
    args = parse_args()

    # Safety: never overwrite existing model v1
    v1_path = os.path.join(os.path.dirname(args.output_model),
                           'site14_site5_combined_model.pkl')
    if os.path.abspath(args.output_model) == os.path.abspath(v1_path):
        print("ERROR: output-model path would overwrite site14_site5_combined_model.pkl! "
              "Choose a different name (e.g. crispr_model_v2.pkl).", file=sys.stderr)
        sys.exit(1)

    if os.path.basename(args.output_model) == 'site14_site5_combined_model.pkl':
        print("ERROR: output filename must not be site14_site5_combined_model.pkl",
              file=sys.stderr)
        sys.exit(1)

    print("Loading training data: {}".format(args.training_data), flush=True)
    if args.neg_pos_ratio is not None:
        print(f"  Streaming load with neg:pos ratio {args.neg_pos_ratio}:1 ...", flush=True)
    df = load_data_balanced(args.training_data, neg_pos_ratio=args.neg_pos_ratio)
    n_pos = df['label'].sum()
    n_neg = (df['label'] == 0).sum()
    print(f"  {len(df):,} rows  ({n_pos:,} positive, {n_neg:,} negative)", flush=True)
    print(f"  {df['individual_id'].nunique()} individuals, {df['dataset'].nunique()} datasets",
          flush=True)

    # Validate features
    missing = [f for f in FEATURE_NAMES if f not in df.columns]
    if missing:
        print(f"ERROR: missing feature columns: {missing}", file=sys.stderr)
        sys.exit(1)

    df = df.dropna(subset=FEATURE_NAMES + ['label'])
    df['label'] = df['label'].astype(int)

    # Data quality checks
    print("\nLabel distribution by dataset:", flush=True)
    grp = df.groupby('dataset')['label'].value_counts().unstack(fill_value=0)
    print(grp.to_string(), flush=True)

    print("\nLabel distribution by individual_id:", flush=True)
    grp2 = df.groupby('individual_id')['label'].value_counts().unstack(fill_value=0)
    print(grp2.to_string(), flush=True)

    # --- Cross-validation ---
    cv_results = None
    if not args.no_cv:
        if df[args.cv_by].nunique() < 2:
            print(f"WARNING: only 1 unique {args.cv_by} — skipping CV", flush=True)
        else:
            cv_results = run_cv(df, args.cv_by, args.threshold,
                                neg_pos_ratio=args.neg_pos_ratio)

    # --- Train final model on all data ---
    print("\nTraining final model on all data ...", flush=True)
    X = df[FEATURE_NAMES].values
    y = df['label'].values

    n_pos = y.sum()
    n_neg = len(y) - n_pos
    spw   = (n_neg / max(n_pos, 1)) if args.neg_pos_ratio is None else 1.0
    print(f"  scale_pos_weight = {spw:.2f}", flush=True)

    final_model = _make_model(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        scale_pos_weight=spw,
    )
    final_model.fit(X, y)
    print("  Done.", flush=True)

    print_feature_importance(final_model)

    # In-sample performance (informational)
    proba_all = final_model.predict_proba(X)[:, 1]
    pred_all  = (proba_all >= args.threshold).astype(int)
    print(f"\nIn-sample AUC: {roc_auc_score(y, proba_all):.4f}", flush=True)
    print(classification_report(y, pred_all, target_names=['background', 'CRISPR'],
                                 zero_division=0), flush=True)

    # --- Save model ---
    os.makedirs(os.path.dirname(os.path.abspath(args.output_model)), exist_ok=True)
    joblib.dump(final_model, args.output_model)
    print(f"\nModel saved: {args.output_model}", flush=True)

    # Verify v1 is still intact
    if os.path.exists(v1_path):
        print(f"Verified: {v1_path} still exists ({os.path.getsize(v1_path):,} bytes)",
              flush=True)

    # Save CV summary if available
    if cv_results is not None:
        summary_path = args.output_model.replace('.pkl', '_cv_summary.csv')
        cv_results['fold_results'].to_csv(summary_path, index=False)
        print(f"CV summary saved: {summary_path}", flush=True)
        print(f"\nFinal cross-validation AUC: {cv_results['overall_auc']:.4f}", flush=True)

    print("\nDone.", flush=True)


if __name__ == '__main__':
    main()
