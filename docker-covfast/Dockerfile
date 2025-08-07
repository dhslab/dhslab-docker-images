# Dockerfile for the covfast tool
# This builds the tool from source using the provided Makefile.

# Use Ubuntu 20.04 (Focal Fossa) as the base image
FROM ubuntu:focal

# Set a working directory inside the container
WORKDIR /tmp/

# Install build dependencies for htslib and covfast
# - build-essential: provides gcc, make, etc.
# - wget: to download htslib
# - zlib1g-dev, libbz2-dev, liblzma-dev, libcurl4-openssl-dev: required by htslib
# - libssl-dev: required for libcrypto dependency
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/dhslab/covfast/archive/refs/tags/v0.1.tar.gz && \
    tar -xvzf v0.1.tar.gz && \
    cd covfast-0.1 && \
    make all && make install && make clean
