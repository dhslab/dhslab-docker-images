# Use the Rocker R image
FROM ghcr.io/dhslab/docker-rbioconductor:3.19

ARG MAMBA_VERSION="24.3.0-0"

ENV CONDA_DIR /opt/conda
COPY env/*.yml /tmp/
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/download/${MAMBA_VERSION}/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -f -b -p /opt/conda && \
    rm -f Miniforge3-$(uname)-$(uname -m).sh && \
    $CONDA_DIR/bin/mamba env update -n base -f /tmp/basic-tools.yml && \
    $CONDA_DIR/bin/conda clean -y --all

# Entry point for the container
CMD ["bash"]