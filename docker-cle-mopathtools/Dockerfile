FROM python:3.11-slim

LABEL maintainer="David Spencer <dspencer@wustl.edu>"
LABEL description="Bioinformatics tools with Python 3.11 and pip-based installs"

ARG DEBIAN_FRONTEND=noninteractive
ARG MOPATHTOOLS_REF=dev
ENV TZ=America/Chicago

WORKDIR /app

# System tools and bioinformatics command-line packages.
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    wget \
    gawk \
    tzdata \
    openssh-client \
    grep \
    samtools \
    bcftools \
    bedtools \
    tabix \
    build-essential \
    && ln -sf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && dpkg-reconfigure --frontend noninteractive tzdata \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies with pip (no conda/mamba).
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install mopathtools package from GitHub.
RUN python -m pip install --no-cache-dir \
    "git+https://github.com/dhslab/cle-mopathtools.git@${MOPATHTOOLS_REF}"

# Extra runtime deps used by shipped scripts but not declared in mopathtools package metadata.
RUN python -m pip install --no-cache-dir \
    pyranges \
    biopython

# Make workflow scripts available on PATH.
RUN set -eux; \
    base_url="https://raw.githubusercontent.com/dhslab/cle-mopathtools/${MOPATHTOOLS_REF}/bin"; \
    for f in \
      filter_cnvs.py \
      filter_and_collect_svs.py \
      create_fastq_list.py \
      prepare_dragen_demux.py \
      parse_input_samplesheet.py \
      merge_and_filter_smallvariants.py \
      filter_sv.py \
      collect_svs_v3.py; do \
      curl -fsSL "${base_url}/${f}" -o "/usr/local/bin/${f}"; \
    done
RUN chmod +x /usr/local/bin/*.py

CMD ["python", "-c", "import mopathtools; print(mopathtools.__version__)"]
