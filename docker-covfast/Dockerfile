# Dockerfile for the covfast tool
# This builds the tool from source using the provided Makefile.

# Use Ubuntu 20.04 (Focal Fossa) as the base image
FROM ubuntu:focal

# Set a working directory inside the container
WORKDIR /opt/

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

# Copy the source code and Makefile into the container
COPY covfast.c .
COPY Makefile .

# Build and install the application using the Makefile's 'install' target
RUN make install && make clean
