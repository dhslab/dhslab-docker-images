# docker-dafqc : self-contained image for nf-dafseq step 2 (the DAF-QC-SMK wrap).
# Bundles pixi + a pinned DAF-QC-SMK checkout + Snakemake 8.21, with the workflow's two conda
# envs (cmd.yaml, python.yaml) PRE-BUILT into a fixed prefix so the step needs no network or
# env-solving at runtime. Published by CI to ghcr.io/dhslab/docker-dafqc:latest.
FROM debian:bookworm-slim

ARG DAFQC_COMMIT=43184be6a44f436d4a32be310fbb43df990a82b2
ARG PIXI_VERSION=0.56.0

LABEL version="1.0"
LABEL description="nf-dafseq step 2: pixi + DAF-QC-SMK (pinned) + prebuilt Snakemake conda envs"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates curl bash procps bc && \
    rm -rf /var/lib/apt/lists/*

# --- pixi (pinned) ---
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v${PIXI_VERSION} PIXI_HOME=/opt/pixi bash
ENV PATH="/opt/pixi/bin:${PATH}"

# --- DAF-QC-SMK pinned at a known-good commit ---
RUN git clone https://github.com/StergachisLab/DAF-QC-SMK.git /opt/DAF-QC-SMK && \
    git -C /opt/DAF-QC-SMK checkout ${DAFQC_COMMIT}

# --- pixi env: snakemake 8.21 + conda (from the repo's pixi.lock) ---
RUN pixi install --manifest-path /opt/DAF-QC-SMK/pixi.toml

# --- pre-build the snakemake conda envs into a fixed, image-local prefix ---
# DAF-QC-SMK builds its DAG purely from the manifest TSV (it does not read the BAM at parse
# time), so a throwaway 1-row manifest + placeholder inputs is enough for Snakemake to enumerate
# every rule's conda env and create it, then exit (--conda-create-envs-only) without running a
# job. Snakemake 8 keys envs by file content, so the prefix matches at runtime -> no rebuild.
ENV SNAKEMAKE_CONDA_PREFIX=/opt/snakemake-conda-envs
RUN mkdir -p /opt/snakemake-conda-envs /opt/prebake && cd /opt/prebake && \
    printf '>chr1\n%s\n' "$(head -c 2000 /dev/zero | tr '\0' 'A')" > ref.fa && \
    : > reads.bam && \
    printf 'sample\tfile\tregs\n'                          >  config.tbl && \
    printf 'prebake\t/opt/prebake/reads.bam\tchr1:1-100\n' >> config.tbl && \
    printf 'ref: /opt/prebake/ref.fa\nmanifest: /opt/prebake/config.tbl\nplatform: ont\nis_fastq: False\n' > config.yaml && \
    pixi run --manifest-path /opt/DAF-QC-SMK/pixi.toml \
        snakemake --configfile config.yaml \
        --software-deployment-method conda \
        --conda-prefix /opt/snakemake-conda-envs \
        --conda-create-envs-only --cores 1 && \
    rm -rf /opt/prebake

# Defaults nf-dafseq points at (overridable by its washu profile back to cluster paths).
ENV DAFQC_REPO=/opt/DAF-QC-SMK
