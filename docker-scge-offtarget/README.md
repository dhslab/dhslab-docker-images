# docker-scge-offtarget

Lightweight image for the **SCGE off-target analysis** step — a pinned Python ML/genomics stack
plus the analysis scripts shipped in [`bin/`](bin/). Unlike [`docker-scge`](../docker-scge) (a
heavy image that builds HTSlib/samtools/bcftools from source and carries samtools/bcftools/
bedtools/nextflow/cyvcf2/biotite/xgboost/nf-core for the *default* SCGE pipeline), this image is a
slim `python:3.11-slim` build scoped to just the off-target ML modules, with dependency versions
locked for reproducibility.

The 8 off-target scripts (`find_edited_reads`, `worklist_from_vcf`, `pon_filter`, `score`,
`hotspot_to_table`, `join_training_table`, `recall_vs_vaf`, `reconcile_offtarget_report` + the
`features`/`pileup_snapshot` helpers) import: pandas, numpy, pysam, scikit-learn, joblib,
matplotlib, scipy, edlib, pyranges.

Pinned stack:

| Package | Version | Notes |
|---|---|---|
| pysam | 0.24.0 | ships cp311 wheels |
| pandas | 2.3.3 | |
| numpy | 2.3.5 | |
| scipy | (unpinned) | imported directly; pin to base-image version once harvested |
| scikit-learn | 1.8.0 | **must match the newest model pickle — see below** |
| joblib | 1.5.2 | |
| matplotlib | 3.10.8 | |
| pyranges | `<1` (0.1.x) | v0 `.df` API required by `find_edited_reads.py`; 1.x removed it |
| edlib | (unpinned) | imported by `find_edited_reads.py`; pin once harvested |
| openpyxl | 3.1.5 | not imported by any off-target script; kept for parity |

Base: `python:3.11-slim` (matches the proven `docker-baseimage` Python) + `build-essential`, so
`pyranges<1`'s Cython deps (`ncls`, `sorted_nearest`) and `edlib` compile if no cp311 wheel exists.

Published by CI to `ghcr.io/dhslab/docker-scge-offtarget` (tags `latest` and a `YYMMDD` date
stamp — see [tagging](#tagging)).

## scikit-learn / model-pickle compatibility

The shipped models were pickled under **different** sklearn versions:

| Model | Used by | Pickled under |
|---|---|---|
| `wgs_shape_model.pkl` | `score.py`, `worklist_from_vcf.py` | sklearn 1.8.0 |
| `site14_site5_combined_model.pkl` | `find_edited_reads.py` (ECS_INDELS, default `--crispr-model`) | sklearn 1.6.1 |

`scikit-learn==1.8.0` loads its own model exactly and **forward**-loads the 1.6.1 model (emits
`InconsistentVersionWarning`, usually fine). Loading the 1.8.0 pickle under an older runtime is the
backward direction and can *raise* (HistGradientBoosting internals changed) — so do not lower this
pin. **Verify both models load** in the build/CI; if the 1.6.1 model errors, re-pickle it under
1.8.0.

## Layout

```
docker-scge-offtarget/
├── Dockerfile
├── README.md
└── bin/            # off-target analysis scripts, placed on PATH at /opt/scge-offtarget/bin
```

## Usage

```bash
docker run --rm -v "$PWD":/data ghcr.io/dhslab/docker-scge-offtarget:latest \
    <script>.py [args...]
```

## Wiring into nf-scge

This image covers **only** the off-target ML scripts — it is *not* a drop-in for
`docker-scge:latest`. The 8 OFFTARGET modules currently declare `container
'ghcr.io/dhslab/docker-scge:latest'`; they must be repointed to
`ghcr.io/dhslab/docker-scge-offtarget:latest`. That change lives in the **nf-scge pipeline repo**,
not here. Do **not** retag this as `docker-scge:latest` — that would break the default SCGE
pipeline and the VEP/report/transgene steps, which need the heavy image's tools.

## Tagging

CI (`.github/workflows/push_docker.yaml`) builds any directory whose `Dockerfile` changed and tags
the image `:latest` and `:YYMMDD` — it does **not** read `LABEL version`. To publish an explicit
`:1.0` tag, add a raw/semver tag to the workflow's `docker/metadata-action` step:

```yaml
tags: |
  type=raw,value=latest,enable=true
  type=raw,value={{date 'YYMMDD'}},enable=true
  type=raw,value=1.0,enable=true
```
