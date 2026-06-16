# docker-dafseq

Image with the tools used by the **nf-dafseq** DAF-seq pipeline (alignment, combinatorial
dedup, and SNP/deletion phasing). Built on `dhspence/docker-baseimage`.

Adds, on top of the base image:

- **minimap2** (pinned to v2.30) for `map-ont` alignment
- **seaborn** for the dedup QC plots (the only plotting dep not already in the base)

The base image already supplies samtools, pysam, pandas, numpy, scipy, scikit-learn, matplotlib,
and the UCSC `bedGraphToBigWig` used to make the phased bigWigs.

Published by CI to `ghcr.io/dhslab/docker-dafseq:latest`.

Covers nf-dafseq steps 1 (`MINIMAP_ALIGN`), 3 (`MARK_DUPLICATES`), and 4 (`PHASE_READS`).
Step 2 (the wrapped DAF-QC-SMK Snakemake) ships as its own image — it bundles pixi + Snakemake
and does not fit this base-image pattern.

## Build

```bash
docker build -t ghcr.io/dhslab/docker-dafseq:latest docker-dafseq
```
