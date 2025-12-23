FROM rocker/r-ver:4.4.0

LABEL description="Container with R and basic packages"

RUN apt-get update && apt-get install -y --no-install-recommends \
    locales \
    libjpeg-turbo8 \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV TERM=xterm

    ## Add a default CRAN mirror
RUN Rscript -e 'install.packages(c("tidyverse","ggplot2","reshape2","cowplot","dplyr","tidyr","qpdf","jsonlite","optparse"))'

# Install bioconductor packages
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' && \
    Rscript -e 'BiocManager::install(c("BiocGenerics","IRanges","S4Vectors","GenomicRanges","XVector","rtracklayer"), ask=FALSE)'

