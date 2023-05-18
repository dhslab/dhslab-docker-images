FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Image with conda and mustache"

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget build-essential libcurl4-openssl-dev && \
    apt-get clean

# Install Conda and modbam2bed
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh && \
    /bin/bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Miniconda3-py39_4.12.0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/pip3 install mustache-hic && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# docker build --no-cache -t mdivr/mustache:v1 .