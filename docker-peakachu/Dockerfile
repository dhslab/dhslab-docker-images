FROM ubuntu:focal-20221130

LABEL version="1.1"
LABEL description="Image with conda and peakachu"

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget build-essential libcurl4-openssl-dev && \
    apt-get clean

# Install conda (mamba version)
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh" && \
    bash Mambaforge-24.3.0-0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Mambaforge-24.3.0-0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/conda config --add channels defaults && \
    $CONDA_DIR/bin/conda config --add channels conda-forge && \
    $CONDA_DIR/bin/conda config --add channels bioconda && \
    $CONDA_DIR/bin/mamba create -n peakachu cooler numba scikit-learn=1.1.2 joblib=1.1.0 && \
    $CONDA_DIR/bin/mamba activate peakachu && \
    $CONDA_DIR/bin/pip3 install -U peakachu && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/


ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

