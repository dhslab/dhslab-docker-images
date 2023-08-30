FROM rocker/ml:4.1.1-cuda11.1

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

ENV CTAN_REPO=https://www.texlive.info/tlnet-archive/2021/10/31/tlnet
ENV QUARTO_VERSION=1.0.36

RUN /rocker_scripts/install_verse.sh
RUN /rocker_scripts/install_quarto.sh
RUN /rocker_scripts/install_geospatial.sh
