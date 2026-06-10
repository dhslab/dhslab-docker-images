# docker-dafqc

Self-contained image for **nf-dafseq** step 2 — the QC + strand-decoration step that wraps the
published [DAF-QC-SMK](https://github.com/StergachisLab/DAF-QC-SMK) Snakemake workflow.

Bundles:

- **pixi** (v0.56.0) + **Snakemake 8.21** (from DAF-QC-SMK's `pixi.lock`)
- a **pinned DAF-QC-SMK** checkout at `/opt/DAF-QC-SMK` (commit `43184be`; bump `ARG DAFQC_COMMIT`)
- the workflow's two conda envs (`workflow/envs/cmd.yaml`, `python.yaml`) **pre-built** into
  `/opt/snakemake-conda-envs`, so the step does no env-solving or network access at runtime

Published by CI to `ghcr.io/dhslab/docker-dafqc:latest`. Companion to
[`docker-dafseq`](../docker-dafseq) (steps 1/3/4); the two cover the whole nf-dafseq pipeline
under a container runtime.

## How the env prebake works

`snakemake --conda-create-envs-only` needs to construct the DAG to know which conda envs to
create. DAF-QC-SMK builds its DAG from the manifest TSV alone (it does not read the input BAM at
parse time), so the build feeds it a throwaway 1-row manifest with placeholder `ref.fa`/`reads.bam`
just to enumerate the rules, creates every env, then exits without running a job. Snakemake 8 keys
conda envs by file content, so the prefix baked here is reused as-is at runtime.

If a future DAF-QC-SMK revision needs valid inputs to build the DAG, swap the prebake `RUN` to use
the repo's own test data instead (`pixi run --manifest-path /opt/DAF-QC-SMK/pixi.toml test-data`,
then run `--conda-create-envs-only` against `test.yaml`).

## Note

This is a heavy build (pixi + Snakemake + two conda envs incl. samtools/minimap2/pbmarkdup/
rustybam and pysam/numpy/pandas/matplotlib/pyabpoa/panel) — expect a multi-GB image and a
several-minute CI build.
