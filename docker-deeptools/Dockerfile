FROM ubuntu:jammy-20240530

ARG MAMBA_VERSION="24.3.0-0"

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget build-essential curl


# Install conda
ENV CONDA_DIR /opt/conda
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/download/${MAMBA_VERSION}/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -f -b -p /opt/conda && \
    rm -f Miniforge3-$(uname)-$(uname -m).sh && \
    $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda bcftools==1.17 htslib==1.17 bedtools=2.31 pandas==2.0.2 pybedtools==0.9.0 openpyxl==3.1.2 samtools==1.18 scipy==1.11.4 deeptools==3.5.5 && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin










# FROM ubuntu:jammy-20240530

# LABEL version="1.0"
# LABEL description="Image with basic bioinformatics tools"

# # Install essential packages
# RUN apt-get update -y && \
#     apt-get install -y wget less

# # Install conda (mamba version) with Jupyter
# ENV CONDA_DIR /opt/conda
# RUN cd /tmp && mkdir -p $CONDA_DIR && \
#     wget "https://github.com/conda-forge/miniforge/releases/download/23.1.0-1/Mambaforge-23.1.0-1-Linux-x86_64.sh" && \
#     bash Mambaforge-23.1.0-1-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
#     rm -f Mambaforge-23.1.0-1-Linux-x86_64.sh && \
#     $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda bcftools==1.17 htslib==1.17 bedtools=2.31 pandas==2.0.2 pybedtools==0.9.0 openpyxl==3.1.2 samtools==1.18 scipy==1.11.4 && \
#     $CONDA_DIR/bin/conda clean -y --all && \
#     ln -s /opt/conda/bin/* /usr/local/bin/

# ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# # docker build --no-cache -t mdivr/basic-tools:v2 .