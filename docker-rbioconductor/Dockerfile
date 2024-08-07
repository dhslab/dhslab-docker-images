# Use the Rocker R image
FROM rocker/r-ver:4.4.0

# Copy installBioc.r script
RUN cp /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r

# Install system dependencies
RUN apt-get update && apt-get -y install \
    less \
    curl \
    gzip \
    lsof \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2 \
    libxml2-dev \
    libxt6 \
    libglpk-dev \
    libmysqlclient21 \
    libcurl4-openssl-dev && \
    apt-get clean

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager'); BiocManager::install(version = '3.19')"

# Install required R and Bioconductor packages
RUN R -e "install.packages(c('argparse', 'data.table', 'dplyr', 'tidyverse', 'vcfR'))"
RUN R -e "BiocManager::install(c('HDF5Array', 'bsseq', 'DSS', 'GenomicRanges', 'IRanges', 'methylKit'))"

ARG MAMBA_VERSION="24.3.0-0"

ENV CONDA_DIR /opt/conda
COPY env/*.yml /tmp/
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/download/${MAMBA_VERSION}/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -f -b -p /opt/conda && \
    rm -f Miniforge3-$(uname)-$(uname -m).sh && \
    $CONDA_DIR/bin/mamba install -y -f /tmp/basic-tools.yml && \
    $CONDA_DIR/bin/conda clean -y --all

# Entry point for the container
CMD ["bash"]