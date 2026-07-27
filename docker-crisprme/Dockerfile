# ghcr.io/dhslab/docker-crisprme
#
# Thin wrapper over the upstream CRISPRme image
FROM pinellolab/crisprme:latest

USER root

RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
