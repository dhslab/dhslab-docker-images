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
    libcurl4-openssl-dev \
    python3 \
    python3-pip && \
    apt-get clean

# Install basic  R packages
RUN R -e "install.packages(c('argparse', 'data.table', 'dplyr', 'tidyverse', 'tidyr', 'stringr', 'purrr', 'readr', 'vcfR'))"

# Install Bioconductor and packages
RUN R -e "install.packages('BiocManager'); BiocManager::install(version = '3.19')"
RUN R -e "BiocManager::install(c('HDF5Array', 'bsseq', 'DSS', 'GenomicRanges', 'IRanges', 'methylKit'))"

# Entry point for the container
CMD ["bash"]