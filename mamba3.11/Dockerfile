# Start from Ubuntu 22.04 LTS for a modern, stable base
FROM ubuntu:22.04

# Add labels for metadata
LABEL maintainer="David Spencer <dspencer@wustl.edu>"
LABEL description="Base Ubuntu 22.04 Image"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# --- 1. Install Basic System Dependencies ---
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    wget \
    gawk \
    less \
    tzdata \
    openssh-client \
    unzip \
    zip \
    bc \
    grep && \
    # Set timezone
    ln -sf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    # Clean up to reduce image size
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR /opt/conda
WORKDIR /tmp
RUN mkdir ${CONDA_DIR} && \
        wget --quiet "https://github.com/conda-forge/miniforge/releases/download/25.3.1-0/Miniforge3-25.3.1-0-Linux-x86_64.sh" -O /tmp/mamba.sh

RUN /bin/bash /tmp/mamba.sh -f -b -p $CONDA_DIR && \
    rm /tmp/mamba.sh

# Add /opt/conda/bin to PATH
ENV PATH=/opt/conda/bin:$PATH

ARG PYTHON_VERSION=3.11

RUN mamba install -y -c conda-forge -c bioconda \
    "python=${PYTHON_VERSION}" && \
    # Clean up conda cache to save space
    mamba clean -y --all