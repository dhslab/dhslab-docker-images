FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Basic image with conda and modbam2bed"

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget && \
    apt-get clean

# Install Conda and modbam2bed
RUN CONDA_DIR=/usr/local && \
    cd /tmp && mkdir -p $CONDA_DIR && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh && \
    /bin/bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Miniconda3-py39_4.12.0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/conda install  -y -c bioconda -c conda-forge -c epi2melabs modbam2bed=0.6.3 && \
    $CONDA_DIR/bin/conda clean -y --all

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin