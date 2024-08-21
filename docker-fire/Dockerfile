FROM ubuntu:jammy-20240530

ARG MAMBA_VERSION="24.3.0-0"

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget rsync git gnutls-bin build-essential curl


# Install ucsc
RUN mkdir  -p /opt/kentutils/bin  && \
    cd /opt/kentutils/bin  && \
    rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./  && \
    ln -s /opt/kentutils/bin/* /usr/local/bin/

# Install conda
ENV CONDA_DIR /opt/conda
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/download/${MAMBA_VERSION}/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -f -b -p /opt/conda && \
    rm -f Miniforge3-$(uname)-$(uname -m).sh && \
    $CONDA_DIR/bin/mamba install -y --channel conda-forge bioconda::snakemake && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin