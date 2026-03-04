FROM python:3.11-slim

LABEL maintainer="David Spencer <dspencer@wustl.edu>"
LABEL description="Mopath bioinformatics tools with Python 3.11"

ARG DEBIAN_FRONTEND=noninteractive
ARG MOPATHTOOLS_REF=dev
ENV TZ=America/Chicago

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
    less \
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

# Extra runtime deps used by shipped scripts but not declared in mopathtools package metadata.
RUN python -m pip install --no-cache-dir \
    pyranges \
    biopython \
    numpy \
    pandas \
    pysam \
    natsort

RUN cp bin/*.py /usr/local/bin/ & chmod +x /usr/local/bin/*.py

