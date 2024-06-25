FROM registry.gitlab.com/quarto-forge/docker/quarto_all

COPY --chown=$MAMBA_USER:$MAMBA_USER polyglot.env.yaml /tmp/env.yaml

COPY --chown=$MAMBA_USER:$MAMBA_USER CoveragePlots.qmd /scripts/CoveragePlots.qmd

RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba install -y -n base -c conda-forge plotnine \
    && micromamba clean --all --yes \
    && rm -rf /opt/conda/conda-meta /tmp/env.yaml

# Set permissions explicitly
RUN ["chmod", "777", "/scripts/CoveragePlots.qmd"]

