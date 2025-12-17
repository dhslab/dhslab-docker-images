FROM mambaorg/micromamba:1.5.8

# Set the environment variable for non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive

# Copy the environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER scge_report_env.yaml /tmp/env.yaml

# Install dependencies into the base environment
# We install into 'base' to avoid needing 'conda activate' logic in Nextflow
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Add the micromamba bin directory to PATH explicitly (ensure quarto is found)
ENV PATH="/opt/conda/bin:$PATH"

# Verify installation
RUN quarto check install

# Set the working directory
WORKDIR /data

