#!/usr/bin/env python3
"""
evaluate_aavs1_v3.py

Held-out evaluation of the CRISPR v3 model on AAVS1 curated data.

Strategy:
  - Stratified 80/20 train/test split (same individual, random seed=42)
  - Retrain XGBoost with v3 hyperparameters on 80%
  - Evaluate on held-out 20%
  - Produce ROC curve, confusion matrix, feature importance, and feature
    distribution plots suitable for PI presentation

Outputs (all in logs/):
  aavs1_v3_eval_roc.png
  aavs1_v3_eval_cm.png
  aavs1_v3_eval_features.png
  aavs1_v3_eval_dist.png
  aavs1_v3_eval_summary.txt
"""

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve,
    confusion_matrix, classification_report,
    average_precision_score,
)
from xgboost import XGBClassifier

PARQUET = os.path.join(REPO, 'training_data', 'aavs1_tptn_only.parquet')
LOG_DIR = os.path.join(REPO, 'logs')
os.makedirs(LOG_DIR, exist_ok=True)

FEATURE_NAMES = [
    'read_pair_gap', 'read_softclip',
    'read_del_vs_control', 'read_ins_vs_control', 'read_mismatch_vs_control',
    'deletion_exclusive_to_edited',
    'is_at_any_target_site', 'distance_to_closest_pam',
    'total_indel_size', 'indel_size_category', 'indel_complexity_score',
    'softclip_dist_to_PAM', 'is_in_homopolymer', 'dist_to_microsatellite',
]

FEATURE_LABELS = {
    'read_pair_gap':              'Read pair gap',
    'read_softclip':              'Read softclip (bp)',
    'read_del_vs_control':        'Del vs control',
    'read_ins_vs_control':        'Ins vs control',
    'read_mismatch_vs_control':   'Mismatch vs control',
    'deletion_exclusive_to_edited': 'Del exclusive to edited',
    'is_at_any_target_site':      'At any target site',
    'distance_to_closest_pam':    'Distance to closest PAM',
    'total_indel_size':           'Total indel size (bp)',
    'indel_size_category':        'Indel size category',
    'indel_complexity_score':     'Indel complexity score',
    'softclip_dist_to_PAM':       'Softclip dist to PAM',
    'is_in_homopolymer':          'In homopolymer',
    'dist_to_microsatellite':     'Dist to microsatellite',
}

# Features displayed in the distribution panel (most biologically interesting)
DIST_FEATURES = [
    'read_mismatch_vs_control', 'read_softclip',
    'read_ins_vs_control', 'read_del_vs_control',
    'indel_complexity_score', 'softclip_dist_to_PAM',
    'dist_to_microsatellite', 'total_indel_size',
]


def main():
    print("=" * 60, flush=True)
    print("CRISPR v3 Model — Held-Out Evaluation (AAVS1)", flush=True)
    print("=" * 60, flush=True)

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------
    print(f"\nLoading {PARQUET} ...", flush=True)
    df = pd.read_parquet(PARQUET)
    df = df.dropna(subset=FEATURE_NAMES + ['label'])
    df['label'] = df['label'].astype(int)
    n_tp = df['label'].sum()
    n_tn = (df['label'] == 0).sum()
    print(f"  {len(df):,} rows  ({n_tp:,} TP  /  {n_tn:,} TN)", flush=True)

    # ------------------------------------------------------------------
    # Stratified 80/20 split
    # ------------------------------------------------------------------
    X = df[FEATURE_NAMES].values
    y = df['label'].values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, random_state=42, stratify=y
    )
    print(f"\nSplit: {len(X_train):,} train  /  {len(X_test):,} test (stratified 80/20)",
          flush=True)
    print(f"  Train: {y_train.sum():,} TP  /  {(y_train == 0).sum():,} TN", flush=True)
    print(f"  Test:  {y_test.sum():,} TP  /  {(y_test == 0).sum():,} TN", flush=True)

    # ------------------------------------------------------------------
    # Train on 80%
    # ------------------------------------------------------------------
    print("\nTraining XGBoost on 80% split ...", flush=True)
    spw = (y_train == 0).sum() / max(y_train.sum(), 1)
    model = XGBClassifier(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        scale_pos_weight=spw,
        eval_metric='logloss', use_label_encoder=False,
        random_state=42, n_jobs=-1, verbosity=0,
    )
    model.fit(X_train, y_train)
    print("  Done.", flush=True)

    # ------------------------------------------------------------------
    # Evaluate on held-out 20%
    # ------------------------------------------------------------------
    proba = model.predict_proba(X_test)[:, 1]
    pred  = (proba >= 0.70).astype(int)

    auc  = roc_auc_score(y_test, proba)
    ap   = average_precision_score(y_test, proba)
    fpr, tpr, _ = roc_curve(y_test, proba)
    prec_curve, rec_curve, _ = precision_recall_curve(y_test, proba)
    cm = confusion_matrix(y_test, pred)

    print(f"\n{'='*60}", flush=True)
    print(f"HELD-OUT TEST RESULTS  (n={len(y_test):,} reads)", flush=True)
    print(f"{'='*60}", flush=True)
    print(f"  ROC-AUC:           {auc:.4f}", flush=True)
    print(f"  Average Precision: {ap:.4f}", flush=True)
    print(flush=True)
    print(classification_report(y_test, pred,
          target_names=['Artifact (TN)', 'True Edit (TP)'], zero_division=0),
          flush=True)
    print(f"Confusion matrix (threshold=0.70):", flush=True)
    print(f"            Pred Artifact  Pred Edit", flush=True)
    print(f"  True Artifact   {cm[0,0]:>8,}  {cm[0,1]:>8,}", flush=True)
    print(f"  True Edit       {cm[1,0]:>8,}  {cm[1,1]:>8,}", flush=True)

    # Feature importance from the held-out model
    imp = pd.Series(model.feature_importances_, index=FEATURE_NAMES).sort_values(ascending=False)
    print("\nFeature importances (held-out model):", flush=True)
    for feat, val in imp.items():
        bar = '#' * int(val * 60)
        print(f"  {feat:<40s} {val:.4f}  {bar}", flush=True)

    # ------------------------------------------------------------------
    # Save text summary
    # ------------------------------------------------------------------
    summary_path = os.path.join(LOG_DIR, 'aavs1_v3_eval_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("CRISPR v3 Held-Out Evaluation — AAVS1\n")
        f.write(f"Dataset: {len(df):,} rows  ({n_tp:,} TP  /  {n_tn:,} TN)\n")
        f.write(f"Split:   {len(X_train):,} train  /  {len(X_test):,} test  (stratified 80/20, seed=42)\n\n")
        f.write(f"ROC-AUC:           {auc:.4f}\n")
        f.write(f"Average Precision: {ap:.4f}\n\n")
        f.write(classification_report(y_test, pred,
                    target_names=['Artifact (TN)', 'True Edit (TP)'], zero_division=0))
        f.write(f"\nConfusion matrix (threshold=0.70):\n")
        f.write(f"            Pred Artifact  Pred Edit\n")
        f.write(f"  True Artifact   {cm[0,0]:>8,}  {cm[0,1]:>8,}\n")
        f.write(f"  True Edit       {cm[1,0]:>8,}  {cm[1,1]:>8,}\n\n")
        f.write("Feature importances (held-out model):\n")
        for feat, val in imp.items():
            f.write(f"  {feat:<40s} {val:.4f}\n")
    print(f"\nText summary: {summary_path}", flush=True)

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        _make_plots(df, fpr, tpr, auc, prec_curve, rec_curve, ap, cm, imp, LOG_DIR)
    except Exception as exc:
        print(f"\nWARNING: Plotting failed ({exc}) — text summary still written.", flush=True)

    print("\nDone.", flush=True)


def _make_plots(df, fpr, tpr, auc, prec_curve, rec_curve, ap, cm, imp, log_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})

    # ------------------------------------------------------------------
    # Figure 1: ROC + Precision-Recall
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        'CRISPR v3 Model — Held-Out Evaluation  (AAVS1, 80/20 stratified split)',
        fontsize=12, fontweight='bold'
    )

    ax = axes[0]
    ax.plot(fpr, tpr, color='steelblue', lw=2.5, label=f'AUC = {auc:.4f}')
    ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', label='Random')
    ax.fill_between(fpr, tpr, alpha=0.08, color='steelblue')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curve')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.02])

    ax = axes[1]
    ax.plot(rec_curve, prec_curve, color='darkorange', lw=2.5,
            label=f'Avg Precision = {ap:.4f}')
    ax.fill_between(rec_curve, prec_curve, alpha=0.08, color='darkorange')
    baseline = df['label'].mean()
    ax.axhline(baseline, color='gray', lw=1, linestyle='--',
               label=f'Baseline ({baseline:.2f})')
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Curve')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.02])

    plt.tight_layout()
    out = os.path.join(log_dir, 'aavs1_v3_eval_roc.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Figure 2: Confusion matrix
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
    plt.colorbar(im, ax=ax)
    classes = ['Artifact\n(TN)', 'True Edit\n(TP)']
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(classes, fontsize=11)
    ax.set_yticklabels(classes, fontsize=11)
    thresh = cm.max() / 2.0
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f'{cm[i, j]:,}',
                    ha='center', va='center', fontsize=14, fontweight='bold',
                    color='white' if cm[i, j] > thresh else 'black')
    ax.set_ylabel('True Label', fontsize=11)
    ax.set_xlabel('Predicted Label (threshold = 0.70)', fontsize=11)
    ax.set_title('Confusion Matrix — Held-Out 20%\nCRISPR v3 Model (AAVS1)', fontsize=12)
    plt.tight_layout()
    out = os.path.join(log_dir, 'aavs1_v3_eval_cm.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Figure 3: Feature importances
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    n = len(imp)
    labels = [FEATURE_LABELS.get(f, f) for f in imp.index]
    vals   = imp.values * 100
    colors = ['steelblue' if v >= 5.0 else '#a8c4d4' for v in vals]

    ax.barh(range(n), vals[::-1], color=colors[::-1])
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels[::-1], fontsize=10)
    ax.set_xlabel('Feature Importance (%)', fontsize=11)
    ax.set_title('Feature Importances — CRISPR v3 (Held-Out Model)', fontsize=12)
    for i, (v, lbl) in enumerate(zip(vals[::-1], labels[::-1])):
        ax.text(v + 0.3, i, f'{v:.1f}%', va='center', fontsize=9)
    ax.grid(True, axis='x', alpha=0.3)
    ax.set_xlim([0, vals.max() * 1.18])
    plt.tight_layout()
    out = os.path.join(log_dir, 'aavs1_v3_eval_features.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Figure 4: Feature distributions by label (2 × 4 grid)
    # ------------------------------------------------------------------
    present = [f for f in DIST_FEATURES if f in df.columns]
    ncols = 4
    nrows = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 8))
    axes = axes.flatten()
    fig.suptitle(
        'Feature Distributions: True Edits vs Artifacts  (AAVS1, n=149,594)\n'
        'Dashed lines = class means  |  Blue = Artifact (TN)  |  Orange = True Edit (TP)',
        fontsize=11, fontweight='bold'
    )

    tn_df = df[df['label'] == 0]
    tp_df = df[df['label'] == 1]

    for i, feat in enumerate(present):
        ax = axes[i]
        tn_vals = tn_df[feat].dropna().values
        tp_vals = tp_df[feat].dropna().values

        # Cap at 99th percentile to avoid sentinel values dominating the axis
        cap = np.percentile(np.concatenate([tn_vals, tp_vals]), 99)
        tn_plot = np.clip(tn_vals, None, cap)
        tp_plot = np.clip(tp_vals, None, cap)

        bins = np.linspace(
            min(tn_plot.min(), tp_plot.min()),
            max(tn_plot.max(), tp_plot.max()),
            60
        )
        ax.hist(tn_plot, bins=bins, alpha=0.55, color='steelblue', density=True,
                label=f'Artifact (n={len(tn_plot):,})')
        ax.hist(tp_plot, bins=bins, alpha=0.55, color='darkorange', density=True,
                label=f'True Edit (n={len(tp_plot):,})')

        tn_mean = float(tn_df[feat].mean())
        tp_mean = float(tp_df[feat].mean())
        ax.axvline(tn_mean, color='steelblue',  lw=2, linestyle='--', alpha=0.9)
        ax.axvline(tp_mean, color='darkorange', lw=2, linestyle='--', alpha=0.9)

        ax.set_title(FEATURE_LABELS.get(feat, feat), fontsize=10, fontweight='bold')
        ax.set_xlabel('Value', fontsize=8)
        ax.set_ylabel('Density', fontsize=8)
        ax.text(0.97, 0.95,
                f'TN μ={tn_mean:.2f}\nTP μ={tp_mean:.2f}',
                transform=ax.transAxes, ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        ax.grid(True, alpha=0.2)

    axes[0].legend(fontsize=8, loc='upper left')

    # Hide any unused cells
    for j in range(len(present), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    out = os.path.join(log_dir, 'aavs1_v3_eval_dist.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)


if __name__ == '__main__':
    main()
