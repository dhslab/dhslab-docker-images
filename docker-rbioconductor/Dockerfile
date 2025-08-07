# Use the Rocker R image
FROM rocker/r-ver:4.4.1

# Copy installBioc.r script
RUN cp /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r

# Install system dependencies
RUN apt-get update && apt-get -y install \
    less \
    curl \
    wget \
    gzip \
    lsof \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2 \
    libxml2-dev \
    libxt6 \
    libglpk-dev \
    libpcre2-dev \
    libdeflate-dev \
    libmysqlclient21 \
    libcurl4-openssl-dev && \
    apt-get clean

# Install basic  R packages
RUN R -e "install.packages(c('argparse', 'data.table', 'dplyr', 'tidyverse', 'tidyr', 'stringr', 'purrr', 'readr', 'vcfR','devtools'))"

# Install Bioconductor and packages
RUN R -e "install.packages('BiocManager'); BiocManager::install(version = '3.19')"
RUN R -e "BiocManager::install(c('HDF5Array', 'bsseq', 'DSS', 'GenomicRanges', 'IRanges', 'methylKit','rtracklayer','DESeq2'))"

# ENV CONDA_DIR /opt/conda
# RUN cd /tmp && mkdir -p $CONDA_DIR && \
#     wget "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh" && \
#     bash Mambaforge-24.3.0-0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
#     rm -f Mambaforge-24.3.0-0-Linux-x86_64.sh && \
#     $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda bcftools==1.20 htslib==1.20 samtools==1.20 bedtools=2.31 pandas==2.0.2 pybedtools==0.9.0 openpyxl==3.1.2 scipy==1.11.4 cython cyvcf2==0.30.28 pysam==0.22.1 statsmodels pyranges==0.0.129 pybigwig biotite && \
#     $CONDA_DIR/bin/conda clean -y --all && \
#     ln -s /opt/conda/bin/* /usr/local/bin/


# Entry point for the container
CMD ["bash"]