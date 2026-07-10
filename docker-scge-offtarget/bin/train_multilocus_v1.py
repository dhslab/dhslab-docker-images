#!/usr/bin/env python3
"""
train_multilocus_v1.py

Train a multi-locus read-level CRISPR edit classifier.

Training data:
  - AAVS1 curated parquet (chr19, 149K reads)
  - v4 parquet (77M reads, 37 samples) minus both holdouts

Holdouts:
  Test A (patient): individual_id == CART_NS0065
                    B2M_2, CTLA4_2, PD1_2, Regnase1_1
                    → tests generalization to a new patient on seen genes

  Test B (gene):    sample_id contains "Regnase1" across all patients
                    NS0027-Regnase1_1/2, NS0065-Regnase1_1, NS0073-Regnase1_2
                    → tests generalization to a completely unseen gene target

  (NS0065-Regnase1_1 appears in both — overlap is intentional)

Training balance:
  All training positives + NEG_RATIO × negatives, sampled proportionally
  per individual_id so no single patient dominates.

Outputs (aavs1_ml_testing/):
  multilocus_v1_summary.txt
  multilocus_v1_roc_{testA,testB}.png
  multilocus_v1_cm_{testA,testB}.png
  multilocus_v1_features.png
  multilocus_v1_by_sample_{testA,testB}.png

Model saved to:
  assets/models/crispr_model_multilocus_v1.pkl
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
from xgboost import XGBClassifier

PARQUET_V4   = os.path.join(REPO, 'training_data', 'crispr_training_v4.parquet')
PARQUET_AAVS1 = os.path.join(REPO, 'training_data', 'aavs1_tptn_only.parquet')
MODEL_OUT    = os.path.join(REPO, 'assets', 'models', 'crispr_model_multilocus_v1.pkl')
OUTPUT_DIR   = os.path.join(REPO, 'aavs1_ml_testing')

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

LOAD_COLS = FEATURE_NAMES + ['label', 'sample_id', 'dataset', 'individual_id']


def is_holdout_A(row_individual_id):
    return row_individual_id == 'CART_NS0065'

def is_holdout_B(row_sample_id):
    return 'Regnase1' in str(row_sample_id)


def stream_and_split(parquet_path, label_override=None):
    """
    Stream parquet, routing each row to train / testA / testB buckets.
    Returns three DataFrames.
    """
    import pyarrow.parquet as pq

    train_chunks, testA_chunks, testB_chunks = [], [], []

    f = pq.ParquetFile(parquet_path)
    for batch in f.iter_batches(batch_size=100_000, columns=LOAD_COLS):
        df = batch.to_pandas()
        df = df.dropna(subset=FEATURE_NAMES + ['label'])
        df['label'] = df['label'].astype(int)
        if label_override is not None:
            df['label'] = label_override

        mask_A = df['individual_id'].apply(is_holdout_A)
        mask_B = df['sample_id'].apply(is_holdout_B)
        mask_train = ~mask_A & ~mask_B

        if mask_A.any():
            testA_chunks.append(df[mask_A])
        if mask_B.any():
            testB_chunks.append(df[mask_B])
        if mask_train.any():
            train_chunks.append(df[mask_train])

    train = pd.concat(train_chunks, ignore_index=True) if train_chunks else pd.DataFrame()
    testA = pd.concat(testA_chunks, ignore_index=True) if testA_chunks else pd.DataFrame()
    testB = pd.concat(testB_chunks, ignore_index=True) if testB_chunks else pd.DataFrame()
    return train, testA, testB


def balance_negatives(df, neg_ratio, seed=42):
    """
    Keep all positives. Sample neg_ratio × n_pos negatives,
    drawn proportionally per individual_id.
    """
    rng = np.random.default_rng(seed)
    pos = df[df['label'] == 1]
    neg = df[df['label'] == 0]

    n_neg_target = len(pos) * neg_ratio
    if len(neg) <= n_neg_target:
        return df

    id_counts = neg['individual_id'].value_counts()
    fracs = (id_counts / len(neg) * n_neg_target).round().astype(int)

    neg_parts = []
    for iid, n in fracs.items():
        pool = neg[neg['individual_id'] == iid]
        n = min(n, len(pool))
        if n > 0:
            idx = rng.choice(len(pool), size=n, replace=False)
            neg_parts.append(pool.iloc[idx])

    neg_sampled = pd.concat(neg_parts, ignore_index=True)
    return pd.concat([pos, neg_sampled], ignore_index=True)


def evaluate(model, df, threshold, label):
    """Run full evaluation on a test DataFrame, return metrics dict + arrays."""
    X = df[FEATURE_NAMES].values
    y = df['label'].values

    proba = model.predict_proba(X)[:, 1]
    pred  = (proba >= threshold).astype(int)

    auc = roc_auc_score(y, proba)
    ap  = average_precision_score(y, proba)
    fpr, tpr, _ = roc_curve(y, proba)
    prec_c, rec_c, _ = precision_recall_curve(y, proba)
    cm = confusion_matrix(y, pred)

    print(f"\n{'='*65}", flush=True)
    print(f"RESULTS — {label}  (n={len(y):,}  TP={y.sum():,}  TN={(y==0).sum():,})",
          flush=True)
    print(f"{'='*65}", flush=True)
    print(f"  ROC-AUC:           {auc:.4f}", flush=True)
    print(f"  Average Precision: {ap:.4f}", flush=True)
    print(flush=True)
    print(classification_report(y, pred,
          target_names=['Artifact (TN)', 'True Edit (TP)'], zero_division=0), flush=True)
    print(f"Confusion matrix (threshold={threshold}):", flush=True)
    print(f"            Pred Artifact  Pred Edit", flush=True)
    print(f"  True Artifact   {cm[0,0]:>8,}  {cm[0,1]:>8,}", flush=True)
    print(f"  True Edit       {cm[1,0]:>8,}  {cm[1,1]:>8,}", flush=True)

    # Per-sample breakdown
    df = df.copy()
    df['_score'] = proba
    df['_pred']  = pred
    sample_stats = []
    for sid in sorted(df['sample_id'].unique()):
        sub = df[df['sample_id'] == sid]
        n_tp = sub['label'].sum()
        n_tn = (sub['label'] == 0).sum()
        if n_tp == 0 or n_tn == 0:
            continue
        sub_auc = roc_auc_score(sub['label'], sub['_score'])
        sub_cm  = confusion_matrix(sub['label'], (sub['_score'] >= threshold).astype(int))
        fp = sub_cm[0, 1]
        tp_ok = sub_cm[1, 1]
        iid = sub['individual_id'].iloc[0]
        print(f"  {sid:<40s}  AUC={sub_auc:.4f}  "
              f"TPR={100*tp_ok/n_tp:.1f}%  FPR={100*fp/n_tn:.2f}%", flush=True)
        sample_stats.append(dict(sample_id=sid, individual_id=iid,
                                 n_tp=n_tp, n_tn=n_tn, roc_auc=sub_auc,
                                 tp_recall=tp_ok/n_tp, fp_rate=fp/n_tn))

    return dict(auc=auc, ap=ap, fpr=fpr, tpr=tpr, prec_c=prec_c, rec_c=rec_c,
                cm=cm, y=y, proba=proba, pred=pred,
                sample_stats=sample_stats, label=label, df=df)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--neg-ratio', type=int, default=5)
    parser.add_argument('--threshold', type=float, default=0.70)
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(MODEL_OUT), exist_ok=True)

    print("=" * 65, flush=True)
    print("CRISPR Multi-Locus Model v1 — Training", flush=True)
    print(f"  Holdout A (patient): CART_NS0065", flush=True)
    print(f"  Holdout B (gene):    all Regnase1 samples", flush=True)
    print(f"  Neg ratio: {args.neg_ratio}:1", flush=True)
    print("=" * 65, flush=True)

    # ------------------------------------------------------------------
    # Stream v4 parquet → train / testA / testB
    # ------------------------------------------------------------------
    print(f"\nStreaming v4 parquet ...", flush=True)
    train_v4, testA, testB = stream_and_split(PARQUET_V4)
    print(f"  v4 train rows:  {len(train_v4):,}  "
          f"(TP={train_v4['label'].sum():,}  TN={(train_v4['label']==0).sum():,})", flush=True)
    print(f"  Test A rows:    {len(testA):,}  "
          f"(TP={testA['label'].sum():,}  TN={(testA['label']==0).sum():,})", flush=True)
    print(f"  Test B rows:    {len(testB):,}  "
          f"(TP={testB['label'].sum():,}  TN={(testB['label']==0).sum():,})", flush=True)

    # ------------------------------------------------------------------
    # Load AAVS1 parquet and tag it
    # ------------------------------------------------------------------
    print(f"\nLoading AAVS1 parquet ...", flush=True)
    aavs1 = pd.read_parquet(PARQUET_AAVS1, columns=LOAD_COLS)
    aavs1 = aavs1.dropna(subset=FEATURE_NAMES + ['label'])
    aavs1['label'] = aavs1['label'].astype(int)
    # Fill in identity columns if missing
    if 'individual_id' not in aavs1.columns or aavs1['individual_id'].isna().all():
        aavs1['individual_id'] = 'AAVS1'
    if 'dataset' not in aavs1.columns or aavs1['dataset'].isna().all():
        aavs1['dataset'] = 'AAVS1_curated'
    print(f"  AAVS1 rows: {len(aavs1):,}  "
          f"(TP={aavs1['label'].sum():,}  TN={(aavs1['label']==0).sum():,})", flush=True)

    # ------------------------------------------------------------------
    # Combine and balance training data
    # ------------------------------------------------------------------
    train_all = pd.concat([train_v4, aavs1], ignore_index=True)
    del train_v4, aavs1  # free memory

    print(f"\nBalancing training data (neg_ratio={args.neg_ratio}) ...", flush=True)
    train_bal = balance_negatives(train_all, neg_ratio=args.neg_ratio)
    del train_all

    train_bal = train_bal.sample(frac=1, random_state=42).reset_index(drop=True)
    n_tp_train = train_bal['label'].sum()
    n_tn_train = (train_bal['label'] == 0).sum()
    print(f"  Balanced training: {len(train_bal):,} rows  "
          f"(TP={n_tp_train:,}  TN={n_tn_train:,})", flush=True)

    iid_dist = train_bal['individual_id'].value_counts()
    print("  Individual distribution in training:", flush=True)
    for iid, cnt in iid_dist.items():
        print(f"    {iid:<30s} {cnt:>8,}", flush=True)

    # ------------------------------------------------------------------
    # Train XGBoost
    # ------------------------------------------------------------------
    print("\nTraining XGBoost ...", flush=True)
    X_train = train_bal[FEATURE_NAMES].values
    y_train = train_bal['label'].values
    del train_bal

    spw = n_tn_train / max(n_tp_train, 1)
    model = XGBClassifier(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        scale_pos_weight=spw,
        eval_metric='logloss', use_label_encoder=False,
        random_state=42, n_jobs=4, verbosity=0,
    )
    model.fit(X_train, y_train)
    del X_train, y_train
    print("  Done.", flush=True)

    # Feature importances
    imp = pd.Series(model.feature_importances_, index=FEATURE_NAMES).sort_values(ascending=False)
    print("\nFeature importances:", flush=True)
    for feat, val in imp.items():
        bar = '#' * int(val * 60)
        print(f"  {feat:<40s} {val:.4f}  {bar}", flush=True)

    # ------------------------------------------------------------------
    # Save model
    # ------------------------------------------------------------------
    artifact = {'model': model, 'threshold': args.threshold, 'features': FEATURE_NAMES}
    joblib.dump(artifact, MODEL_OUT)
    print(f"\nModel saved: {MODEL_OUT}", flush=True)

    # ------------------------------------------------------------------
    # Evaluate on both test sets
    # ------------------------------------------------------------------
    # Balance test sets for evaluation (same 5:1 ratio for tractable scoring)
    testA_bal = balance_negatives(testA, neg_ratio=args.neg_ratio)
    testB_bal = balance_negatives(testB, neg_ratio=args.neg_ratio)

    res_A = evaluate(model, testA_bal, args.threshold, 'Test A — Patient Holdout (CART_NS0065)')
    res_B = evaluate(model, testB_bal, args.threshold, 'Test B — Gene Holdout (Regnase1)')

    # ------------------------------------------------------------------
    # Write text summary
    # ------------------------------------------------------------------
    summary_path = os.path.join(OUTPUT_DIR, 'multilocus_v1_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("CRISPR Multi-Locus Model v1 — Training & Evaluation\n")
        f.write(f"Neg ratio: {args.neg_ratio}:1  |  Threshold: {args.threshold}\n\n")
        f.write(f"Training: {n_tp_train+n_tn_train:,} rows  "
                f"(TP={n_tp_train:,}  TN={n_tn_train:,})\n\n")

        for res, tag in [(res_A, 'Test A — Patient Holdout (CART_NS0065)'),
                         (res_B, 'Test B — Gene Holdout (Regnase1)')]:
            y, pred = res['y'], res['pred']
            cm = res['cm']
            f.write(f"{'='*65}\n{tag}\n{'='*65}\n")
            f.write(f"ROC-AUC:           {res['auc']:.4f}\n")
            f.write(f"Average Precision: {res['ap']:.4f}\n\n")
            f.write(classification_report(y, pred,
                        target_names=['Artifact (TN)', 'True Edit (TP)'], zero_division=0))
            f.write(f"\nConfusion matrix (threshold={args.threshold}):\n")
            f.write(f"  True Artifact   {cm[0,0]:>8,}  {cm[0,1]:>8,}\n")
            f.write(f"  True Edit       {cm[1,0]:>8,}  {cm[1,1]:>8,}\n\n")
            f.write("Per-sample breakdown:\n")
            for row in res['sample_stats']:
                f.write(f"  {row['sample_id']:<40s}  AUC={row['roc_auc']:.4f}  "
                        f"TPR={100*row['tp_recall']:.1f}%  "
                        f"FPR={100*row['fp_rate']:.2f}%\n")
            f.write("\n")

        f.write("Feature importances:\n")
        for feat, val in imp.items():
            f.write(f"  {feat:<40s} {val:.4f}\n")

    print(f"\nSummary: {summary_path}", flush=True)

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        _make_plots(res_A, res_B, imp, args.threshold, OUTPUT_DIR)
    except Exception as exc:
        print(f"\nWARNING: Plotting failed ({exc})", flush=True)
        import traceback; traceback.print_exc()

    print("\nDone.", flush=True)


def _make_plots(res_A, res_B, imp, threshold, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})

    palette = ['steelblue', 'darkorange', 'green', 'purple', 'red']

    # ------------------------------------------------------------------
    # ROC + PR for each test set
    # ------------------------------------------------------------------
    for res, tag in [(res_A, 'testA'), (res_B, 'testB')]:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f'Multi-Locus Model v1 — {res["label"]}',
                     fontsize=11, fontweight='bold')

        ax = axes[0]
        ax.plot(res['fpr'], res['tpr'], color='steelblue', lw=2.5,
                label=f'AUC = {res["auc"]:.4f}')
        ax.plot([0,1],[0,1], color='gray', lw=1, linestyle='--', label='Random')
        ax.axhline(1.0, color='red', linestyle=':', lw=1, alpha=0.5,
                   label='AAVS1 AUC (1.0000)')
        ax.fill_between(res['fpr'], res['tpr'], alpha=0.08, color='steelblue')
        ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curve'); ax.legend(fontsize=9, loc='lower right')
        ax.grid(True, alpha=0.3); ax.set_xlim([0,1]); ax.set_ylim([0,1.02])

        ax = axes[1]
        ax.plot(res['rec_c'], res['prec_c'], color='darkorange', lw=2.5,
                label=f'Avg Precision = {res["ap"]:.4f}')
        ax.fill_between(res['rec_c'], res['prec_c'], alpha=0.08, color='darkorange')
        baseline = res['y'].mean()
        ax.axhline(baseline, color='gray', lw=1, linestyle='--',
                   label=f'Baseline ({baseline:.3f})')
        ax.set_xlabel('Recall'); ax.set_ylabel('Precision')
        ax.set_title('Precision-Recall Curve'); ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3); ax.set_xlim([0,1]); ax.set_ylim([0,1.02])

        plt.tight_layout()
        out = os.path.join(out_dir, f'multilocus_v1_roc_{tag}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
        print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Confusion matrices
    # ------------------------------------------------------------------
    for res, tag in [(res_A, 'testA'), (res_B, 'testB')]:
        fig, ax = plt.subplots(figsize=(6, 5))
        cm = res['cm']
        im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
        plt.colorbar(im, ax=ax)
        classes = ['Artifact\n(TN)', 'True Edit\n(TP)']
        ax.set_xticks([0,1]); ax.set_yticks([0,1])
        ax.set_xticklabels(classes, fontsize=11)
        ax.set_yticklabels(classes, fontsize=11)
        thresh = cm.max() / 2.0
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f'{cm[i,j]:,}', ha='center', va='center',
                        fontsize=13, fontweight='bold',
                        color='white' if cm[i,j] > thresh else 'black')
        ax.set_ylabel('True Label', fontsize=11)
        ax.set_xlabel(f'Predicted (threshold={threshold})', fontsize=11)
        ax.set_title(f'Confusion Matrix\n{res["label"]}', fontsize=10)
        plt.tight_layout()
        out = os.path.join(out_dir, f'multilocus_v1_cm_{tag}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
        print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Feature importances
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    n = len(imp)
    labels = [FEATURE_LABELS.get(f, f) for f in imp.index]
    vals   = imp.values * 100
    colors = ['steelblue' if v >= 5.0 else '#a8c4d4' for v in vals]
    ax.barh(range(n), vals[::-1], color=colors[::-1])
    ax.set_yticks(range(n)); ax.set_yticklabels(labels[::-1], fontsize=10)
    ax.set_xlabel('Feature Importance (%)', fontsize=11)
    ax.set_title('Feature Importances — Multi-Locus Model v1', fontsize=12)
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.3, i, f'{v:.1f}%', va='center', fontsize=9)
    ax.grid(True, axis='x', alpha=0.3); ax.set_xlim([0, vals.max() * 1.18])
    plt.tight_layout()
    out = os.path.join(out_dir, 'multilocus_v1_features.png')
    plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
    print(f"Saved: {out}", flush=True)

    # ------------------------------------------------------------------
    # Per-sample TPR vs FPR scatter — one panel per test set
    # ------------------------------------------------------------------
    for res, tag in [(res_A, 'testA'), (res_B, 'testB')]:
        stats = res['sample_stats']
        if not stats:
            continue
        fig, ax = plt.subplots(figsize=(9, 6))
        iids = sorted(set(r['individual_id'] for r in stats))
        color_map = {iid: palette[i % len(palette)] for i, iid in enumerate(iids)}

        for row in stats:
            color = color_map[row['individual_id']]
            ax.scatter(row['fp_rate'] * 100, row['tp_recall'] * 100,
                       color=color, s=90, alpha=0.85, zorder=3)
            ax.annotate(row['sample_id'].split('-', 1)[-1],
                        (row['fp_rate']*100, row['tp_recall']*100),
                        textcoords='offset points', xytext=(5, 3), fontsize=7, alpha=0.8)

        # Legend by individual
        for iid, color in color_map.items():
            ax.scatter([], [], color=color, s=60, label=iid)
        ax.legend(fontsize=8, loc='lower right')

        ax.set_xlabel('False Positive Rate (%)', fontsize=11)
        ax.set_ylabel('True Positive Rate (%)', fontsize=11)
        ax.set_title(f'Per-Sample Performance — {res["label"]}', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0); ax.set_ylim([0, 105])
        plt.tight_layout()
        out = os.path.join(out_dir, f'multilocus_v1_by_sample_{tag}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight'); plt.close()
        print(f"Saved: {out}", flush=True)


if __name__ == '__main__':
    main()
