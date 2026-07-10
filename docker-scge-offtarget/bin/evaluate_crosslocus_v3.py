#!/usr/bin/env python3
"""
evaluate_crosslocus_v3.py

Cross-locus generalization test for the AAVS1-trained CRISPR v3 model.

Strategy:
  - Load pre-trained model from assets/models/crispr_model_aavs1_v3.pkl
    (XGBoost trained on AAVS1 chr19 reads only)
  - Load v4 parquet (77M reads: CART_wgs, CART_KO_DNA, NS0073_wgs —
    completely different loci: B2M, CTLA4, PD1, Regnase1, KO-DNA genes)
  - Sample: all positives + up to MAX_NEG negatives (stratified by sample_id)
  - Score with AAVS1 model, no retraining
  - Report overall ROC-AUC, PR-AUC, and per-dataset / per-sample breakdown

Outputs (in OUTPUT_DIR):
  crosslocus_v3_roc.png
  crosslocus_v3_cm.png
  crosslocus_v3_features.png
  crosslocus_v3_by_dataset.png
  crosslocus_v3_summary.txt

Usage:
    python bin/evaluate_crosslocus_v3.py [--max-neg N] [--threshold T]
"""

import argparse
import os
import sys
import warnings

warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve,
    confusion_matrix, classification_report,
    average_precision_score,
)

PARQUET_V4    = os.path.join(REPO, 'training_data', 'crispr_training_v4.parquet')
MODEL_PATH    = os.path.join(REPO, 'assets', 'models', 'crispr_model_aavs1_v3.pkl')
OUTPUT_DIR    = os.path.join(REPO, 'aavs1_ml_testing')

FEATURE_NAMES = [
    'read_pair_gap', 'read_softclip',
    'read_del_vs_control', 'read_ins_vs_control', 'read_mismatch_vs_control',
    'deletion_exclusive_to_edited',
    'is_at_any_target_site', 'distance_to_closest_pam',
    'total_indel_size', 'indel_size_category', 'indel_complexity_score',
    'softclip_dist_to_PAM', 'is_in_homopolymer', 'dist_to_microsatellite',
]

FEATURE_LABELS = {
    'read_pair_gap':               'Read pair gap',
    'read_softclip':               'Read softclip (bp)',
    'read_del_vs_control':         'Del vs control',
    'read_ins_vs_control':         'Ins vs control',
    'read_mismatch_vs_control':    'Mismatch vs control',
    'deletion_exclusive_to_edited':'Del exclusive to edited',
    'is_at_any_target_site':       'At any target site',
    'distance_to_closest_pam':     'Distance to closest PAM',
    'total_indel_size':            'Total indel size (bp)',
    'indel_size_category':         'Indel size category',
    'indel_complexity_score':      'Indel complexity score',
    'softclip_dist_to_PAM':        'Softclip dist to PAM',
    'is_in_homopolymer':           'In homopolymer',
    'dist_to_microsatellite':      'Dist to microsatellite',
}


def load_v4_sample(parquet_path, max_neg, seed=42):
    """
    Load all positives + up to max_neg negatives from v4 parquet.
    Negatives are sampled proportionally per sample_id to keep representation balanced.
    """
    import pyarrow.parquet as pq

    rng = np.random.default_rng(seed)
    print(f"Streaming v4 parquet: {parquet_path}", flush=True)
    f = pq.ParquetFile(parquet_path)

    needed_cols = FEATURE_NAMES + ['label', 'sample_id', 'dataset', 'chrom', 'pos']

    pos_chunks = []
    neg_chunks = []

    for batch in f.iter_batches(batch_size=100_000, columns=needed_cols):
        df_b = batch.to_pandas()
        df_b = df_b.dropna(subset=FEATURE_NAMES + ['label'])
        df_b['label'] = df_b['label'].astype(int)
        pos_chunks.append(df_b[df_b['label'] == 1])
        neg_chunks.append(df_b[df_b['label'] == 0])

    pos_df = pd.concat(pos_chunks, ignore_index=True)
    neg_df = pd.concat(neg_chunks, ignore_index=True)

    print(f"  Total positives: {len(pos_df):,}", flush=True)
    print(f"  Total negatives: {len(neg_df):,}", flush=True)

    if len(neg_df) > max_neg:
        # Sample proportionally by sample_id
        sample_counts = neg_df['sample_id'].value_counts()
        fracs = (sample_counts / len(neg_df) * max_neg).round().astype(int)
        neg_parts = []
        for sid, n in fracs.items():
            pool = neg_df[neg_df['sample_id'] == sid]
            n = min(n, len(pool))
            if n > 0:
                idx = rng.choice(len(pool), size=n, replace=False)
                neg_parts.append(pool.iloc[idx])
        neg_sampled = pd.concat(neg_parts, ignore_index=True)
    else:
        neg_sampled = neg_df

    df = pd.concat([pos_df, neg_sampled], ignore_index=True)
    df = df.sample(frac=1, random_state=seed).reset_index(drop=True)
    print(f"  Sampled: {len(pos_df):,} TP + {len(neg_sampled):,} TN = {len(df):,} total",
          flush=True)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-neg', type=int, default=500_000,
                        help='Max negatives to sample from v4 (default: 500000)')
    parser.add_argument('--threshold', type=float, default=0.70,
                        help='Classification threshold (default: 0.70, same as AAVS1 eval)')
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 65, flush=True)
    print("CRISPR v3 Model — Cross-Locus Generalization Test", flush=True)
    print("  Train: AAVS1 (chr19)  →  Test: B2M/CTLA4/PD1/Regnase1/KO-DNA", flush=True)
    print("=" * 65, flush=True)

    # ------------------------------------------------------------------
    # Load model
    # ------------------------------------------------------------------
    print(f"\nLoading model: {MODEL_PATH}", flush=True)
    artifact = joblib.load(MODEL_PATH)
    if isinstance(artifact, dict):
        model = artifact['model']
        saved_threshold = artifact.get('threshold', args.threshold)
        saved_features  = artifact.get('features', FEATURE_NAMES)
        print(f"  Model type: {type(model).__name__}", flush=True)
        print(f"  Saved threshold: {saved_threshold}", flush=True)
        print(f"  Features ({len(saved_features)}): {saved_features}", flush=True)
        threshold = saved_threshold
        feature_names = saved_features
    else:
        model = artifact
        threshold = args.threshold
        feature_names = FEATURE_NAMES
        print(f"  Model type: {type(model).__name__} (bare model)", flush=True)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    df = load_v4_sample(PARQUET_V4, max_neg=args.max_neg)

    X = df[feature_names].values
    y = df['label'].values

    # ------------------------------------------------------------------
    # Score (no retraining — true generalization test)
    # ------------------------------------------------------------------
    print("\nScoring with AAVS1-trained model (no retraining) ...", flush=True)
    proba = model.predict_proba(X)[:, 1]
    pred  = (proba >= threshold).astype(int)
    df['score'] = proba
    df['pred']  = pred

    auc = roc_auc_score(y, proba)
    ap  = average_precision_score(y, proba)
    fpr, tpr, _ = roc_curve(y, proba)
    prec_curve, rec_curve, _ = precision_recall_curve(y, proba)
    cm = confusion_matrix(y, pred)

    # ------------------------------------------------------------------
    # Overall results
    # ------------------------------------------------------------------
    print(f"\n{'='*65}", flush=True)
    print(f"CROSS-LOCUS TEST RESULTS  (n={len(y):,} reads)", flush=True)
    print(f"  TP: {y.sum():,}  |  TN: {(y==0).sum():,}", flush=True)
    print(f"{'='*65}", flush=True)
    print(f"  ROC-AUC:           {auc:.4f}", flush=True)
    print(f"  Average Precision: {ap:.4f}", flush=True)
    print(flush=True)
    print(classification_report(y, pred,
          target_names=['Artifact (TN)', 'True Edit (TP)'], zero_division=0),
          flush=True)
    print(f"Confusion matrix (threshold={threshold}):", flush=True)
    print(f"            Pred Artifact  Pred Edit", flush=True)
    print(f"  True Artifact   {cm[0,0]:>8,}  {cm[0,1]:>8,}", flush=True)
    print(f"  True Edit       {cm[1,0]:>8,}  {cm[1,1]:>8,}", flush=True)

    # ------------------------------------------------------------------
    # Per-dataset breakdown
    # ------------------------------------------------------------------
    print(f"\n{'='*65}", flush=True)
    print("PER-DATASET BREAKDOWN", flush=True)
    print(f"{'='*65}", flush=True)
    dataset_stats = []
    for ds in sorted(df['dataset'].unique()):
        sub = df[df['dataset'] == ds]
        n_tp = sub['label'].sum()
        n_tn = (sub['label'] == 0).sum()
        if n_tp == 0 or n_tn == 0:
            print(f"  {ds}: skipped (no positives or negatives)", flush=True)
            continue
        sub_auc = roc_auc_score(sub['label'], sub['score'])
        sub_ap  = average_precision_score(sub['label'], sub['score'])
        sub_pred = (sub['score'] >= threshold).astype(int)
        sub_cm  = confusion_matrix(sub['label'], sub_pred)
        fp = sub_cm[0, 1]
        fn = sub_cm[1, 0]
        tp_correct = sub_cm[1, 1]
        print(f"\n  [{ds}]  n={len(sub):,}  TP={n_tp:,}  TN={n_tn:,}", flush=True)
        print(f"    ROC-AUC: {sub_auc:.4f}   Avg Precision: {sub_ap:.4f}", flush=True)
        print(f"    TP correctly called: {tp_correct}/{n_tp}  "
              f"({100*tp_correct/n_tp:.1f}%)", flush=True)
        print(f"    FP (artifacts called as edits): {fp}/{n_tn} "
              f"({100*fp/n_tn:.2f}%)", flush=True)
        dataset_stats.append({
            'dataset': ds, 'n': len(sub), 'n_tp': n_tp, 'n_tn': n_tn,
            'roc_auc': sub_auc, 'avg_precision': sub_ap,
            'tp_recall': tp_correct / n_tp,
            'fp_rate': fp / n_tn,
        })

    # ------------------------------------------------------------------
    # Per-sample breakdown
    # ------------------------------------------------------------------
    print(f"\n{'='*65}", flush=True)
    print("PER-SAMPLE BREAKDOWN", flush=True)
    print(f"{'='*65}", flush=True)
    sample_stats = []
    for sid in sorted(df['sample_id'].unique()):
        sub = df[df['sample_id'] == sid]
        n_tp = sub['label'].sum()
        n_tn = (sub['label'] == 0).sum()
        if n_tp == 0 or n_tn == 0:
            continue
        sub_auc = roc_auc_score(sub['label'], sub['score'])
        sub_pred = (sub['score'] >= threshold).astype(int)
        sub_cm  = confusion_matrix(sub['label'], sub_pred)
        fp = sub_cm[0, 1]
        tp_correct = sub_cm[1, 1]
        ds = sub['dataset'].iloc[0]
        print(f"  {sid:<35s} [{ds}]  "
              f"AUC={sub_auc:.4f}  "
              f"TPR={100*tp_correct/n_tp:.1f}%  "
              f"FPR={100*fp/n_tn:.2f}%", flush=True)
        sample_stats.append({
            'sample_id': sid, 'dataset': ds,
            'n_tp': n_tp, 'n_tn': n_tn,
            'roc_auc': sub_auc,
            'tp_recall': tp_correct / n_tp,
            'fp_rate': fp / n_tn,
        })

    # ------------------------------------------------------------------
    # Feature importances (from the loaded model)
    # ------------------------------------------------------------------
    imp = pd.Series(model.feature_importances_, index=feature_names).sort_values(ascending=False)
    print(f"\n{'='*65}", flush=True)
    print("FEATURE IMPORTANCES (AAVS1-trained model)", flush=True)
    print(f"{'='*65}", flush=True)
    for feat, val in imp.items():
        bar = '#' * int(val * 60)
        print(f"  {feat:<40s} {val:.4f}  {bar}", flush=True)

    # ------------------------------------------------------------------
    # Save text summary
    # ------------------------------------------------------------------
    summary_path = os.path.join(OUTPUT_DIR, 'crosslocus_v3_summary.txt')
    with open(summary_path, 'w') as fout:
        fout.write("CRISPR v3 Model — Cross-Locus Generalization Test\n")
        fout.write("Train locus: AAVS1 (chr19)\n")
        fout.write("Test loci:   B2M, CTLA4, PD1, Regnase1, KO-DNA genes\n")
        fout.write(f"Test reads:  {len(df):,}  ({y.sum():,} TP / {(y==0).sum():,} TN)\n\n")
        fout.write(f"ROC-AUC:           {auc:.4f}\n")
        fout.write(f"Average Precision: {ap:.4f}\n\n")
        fout.write(classification_report(y, pred,
                       target_names=['Artifact (TN)', 'True Edit (TP)'], zero_division=0))
        fout.write(f"\nConfusion matrix (threshold={threshold}):\n")
        fout.write(f"            Pred Artifact  Pred Edit\n")
        fout.write(f"  True Artifact   {cm[0,0]:>8,}  {cm[0,1]:>8,}\n")
        fout.write(f"  True Edit       {cm[1,0]:>8,}  {cm[1,1]:>8,}\n\n")
        fout.write("Per-dataset breakdown:\n")
        for row in dataset_stats:
            fout.write(f"  {row['dataset']:<20s}  AUC={row['roc_auc']:.4f}  "
                       f"TPR={100*row['tp_recall']:.1f}%  FPR={100*row['fp_rate']:.2f}%\n")
        fout.write("\nPer-sample breakdown:\n")
        for row in sample_stats:
            fout.write(f"  {row['sample_id']:<35s}  AUC={row['roc_auc']:.4f}  "
                       f"TPR={100*row['tp_recall']:.1f}%  FPR={100*row['fp_rate']:.2f}%\n")
        fout.write("\nFeature importances:\n")
        for feat, val in imp.items():
            fout.write(f"  {feat:<40s} {val:.4f}\n")

    print(f"\nText summary: {summary_path}", flush=True)

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        _make_plots(df, fpr, tpr, auc, prec_curve, rec_curve, ap,
                    cm, imp, threshold, dataset_stats, sample_stats, OUTPUT_DIR)
    except Exception as exc:
        print(f"\nWARNING: Plotting failed ({exc}) — text summary still written.", flush=True)
        import traceback; traceback.print_exc()

    print("\nDone.", flush=True)


def _make_plots(df, fpr, tpr, auc, prec_curve, rec_curve, ap,
                cm, imp, threshold, dataset_stats, sample_stats, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})
    FEATURE_NAMES_LOCAL = list(imp.index)

    # ------------------------------------------------------------------
    # Figure 1: ROC + Precision-Recall
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        'CRISPR v3 — Cross-Locus Generalization\n'
        'Model trained on AAVS1 (chr19), tested on B2M/CTLA4/PD1/Regnase1/KO-DNA',
        fontsize=11, fontweight='bold'
    )

    ax = axes[0]
    ax.plot(fpr, tpr, color='steelblue', lw=2.5, label=f'AUC = {auc:.4f}')
    ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', label='Random')
    ax.fill_between(fpr, tpr, alpha=0.08, color='steelblue')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curve (Cross-Locus)')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1.02])

    ax = axes[1]
    ax.plot(rec_curve, prec_curve, color='darkorange', lw=2.5,
            label=f'Avg Precision = {ap:.4f}')
    ax.fill_between(rec_curve, prec_curve, alpha=0.08, color='darkorange')
    baseline = df['label'].mean()
    ax.axhline(baseline, color='gray', lw=1, linestyle='--',
               label=f'Baseline ({baseline:.4f})')
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Curve (Cross-Locus)')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1.02])

    plt.tight_layout()
    out = os.path.join(out_dir, 'crosslocus_v3_roc.png')
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
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(classes, fontsize=11)
    ax.set_yticklabels(classes, fontsize=11)
    thresh = cm.max() / 2.0
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f'{cm[i,j]:,}', ha='center', va='center',
                    fontsize=14, fontweight='bold',
                    color='white' if cm[i,j] > thresh else 'black')
    ax.set_ylabel('True Label', fontsize=11)
    ax.set_xlabel(f'Predicted Label (threshold = {threshold})', fontsize=11)
    ax.set_title('Confusion Matrix — Cross-Locus\n(AAVS1 model → B2M/CTLA4/PD1/Regnase1/KO-DNA)',
                 fontsize=11)
    plt.tight_layout()
    out = os.path.join(out_dir, 'crosslocus_v3_cm.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Figure 3: Feature importances (model's own weights, unchanged)
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
    ax.set_title('Feature Importances — AAVS1-trained model\n(unchanged; shown for reference)',
                 fontsize=11)
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.3, i, f'{v:.1f}%', va='center', fontsize=9)
    ax.grid(True, axis='x', alpha=0.3)
    ax.set_xlim([0, vals.max() * 1.18])
    plt.tight_layout()
    out = os.path.join(out_dir, 'crosslocus_v3_features.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Figure 4: Per-dataset AUC bar chart + per-sample scatter
    # ------------------------------------------------------------------
    if not dataset_stats and not sample_stats:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Cross-Locus Generalization — Per-Dataset & Per-Sample Performance',
                 fontsize=12, fontweight='bold')

    # Left: per-dataset AUC bars
    ax = axes[0]
    ds_names = [r['dataset'] for r in dataset_stats]
    ds_aucs  = [r['roc_auc'] for r in dataset_stats]
    colors_ds = ['steelblue', 'darkorange', 'green'][:len(ds_names)]
    bars = ax.bar(ds_names, ds_aucs, color=colors_ds, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, color='gray', linestyle='--', lw=1, label='Random')
    ax.axhline(1.0, color='black', linestyle=':', lw=0.8, alpha=0.4)
    # Add AAVS1 reference line
    ax.axhline(1.0, color='red', linestyle='--', lw=1.5, alpha=0.7,
               label='AAVS1 train AUC (1.0000)')
    for bar, val in zip(bars, ds_aucs):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.005,
                f'{val:.4f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax.set_ylim([0, 1.08])
    ax.set_ylabel('ROC-AUC', fontsize=11)
    ax.set_title('ROC-AUC by Dataset', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)

    # Right: per-sample scatter (TPR vs FPR)
    ax = axes[1]
    ds_color_map = {}
    palette = ['steelblue', 'darkorange', 'green', 'purple', 'red']
    for i, row in enumerate(dataset_stats):
        ds_color_map[row['dataset']] = palette[i % len(palette)]

    for row in sample_stats:
        color = ds_color_map.get(row['dataset'], 'gray')
        ax.scatter(row['fp_rate'] * 100, row['tp_recall'] * 100,
                   color=color, s=80, alpha=0.85, zorder=3,
                   label=row['dataset'])

    # Deduplicate legend
    handles, lbls = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, lbls):
        if l not in seen:
            seen[l] = h
    ax.legend(seen.values(), seen.keys(), fontsize=9, loc='lower right')

    ax.axhline(100, color='gray', linestyle=':', lw=0.8, alpha=0.5)
    ax.set_xlabel('False Positive Rate (%)', fontsize=11)
    ax.set_ylabel('True Positive Rate (%)', fontsize=11)
    ax.set_title('TP Rate vs FP Rate by Sample\n(each point = one sample_id)', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim([0, 105])

    plt.tight_layout()
    out = os.path.join(out_dir, 'crosslocus_v3_by_dataset.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}", flush=True)


if __name__ == '__main__':
    main()
