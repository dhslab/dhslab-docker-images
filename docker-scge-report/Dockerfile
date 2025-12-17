FROM mambaorg/micromamba:1.5.8

# Set the environment variable for non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive

# Copy the environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER scge_report_env.yaml /tmp/env.yaml

# Install dependencies into the base environment
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# IMPORTANT: This ARG is required for the "RUN" commands below to see the environment
# It tells the base image's entrypoint to activate the 'base' environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Verify installation
RUN quarto check install

# Set the PATH for the final image so nextflow doesn't need 'conda activate'
ENV PATH="/opt/conda/bin:$PATH"

# Set the working directory
WORKDIR /data
