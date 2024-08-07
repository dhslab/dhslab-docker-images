FROM mdivr/r-ver_4.3.0:v1

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2023.03.1+446
# ENV DEFAULT_USER=rstudio
ARG DEFAULT_USER
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default

COPY scripts/config_R_cuda.sh /rocker_scripts/config_R_cuda.sh 
COPY scripts/default_user.sh /rocker_scripts/default_user.sh 
COPY scripts/init_set_env.sh /rocker_scripts/init_set_env.sh 
COPY scripts/init_userconf.sh /rocker_scripts/init_userconf.sh 
COPY scripts/install_cuda-10.1.sh /rocker_scripts/install_cuda-10.1.sh 
COPY scripts/install_cuda-11.1.sh /rocker_scripts/install_cuda-11.1.sh 
COPY scripts/install_geospatial.sh /rocker_scripts/install_geospatial.sh 
COPY scripts/install_julia.sh /rocker_scripts/install_julia.sh 
COPY scripts/install_jupyter.sh /rocker_scripts/install_jupyter.sh 
COPY scripts/install_nvtop.sh /rocker_scripts/install_nvtop.sh 
COPY scripts/install_pandoc.sh /rocker_scripts/install_pandoc.sh 
COPY scripts/install_pyenv.sh /rocker_scripts/install_pyenv.sh 
COPY scripts/install_python.sh /rocker_scripts/install_python.sh 
COPY scripts/install_quarto.sh /rocker_scripts/install_quarto.sh 
COPY scripts/install_R_ppa.sh /rocker_scripts/install_R_ppa.sh 
COPY scripts/install_R_source.sh /rocker_scripts/install_R_source.sh 
COPY scripts/install_rstudio.sh /rocker_scripts/install_rstudio.sh 
COPY scripts/install_s6init.sh /rocker_scripts/install_s6init.sh 
COPY scripts/install_shiny_server.sh /rocker_scripts/install_shiny_server.sh 
COPY scripts/install_tensorflow.sh /rocker_scripts/install_tensorflow.sh 
COPY scripts/install_texlive.sh /rocker_scripts/install_texlive.sh 
COPY scripts/install_tf1_cuda_10_0.sh /rocker_scripts/install_tf1_cuda_10_0.sh 
COPY scripts/install_tidyverse.sh /rocker_scripts/install_tidyverse.sh 
COPY scripts/install_verse.sh /rocker_scripts/install_verse.sh 
COPY scripts/install_wgrib2.sh /rocker_scripts/install_wgrib2.sh 
COPY scripts/pam-helper.sh /rocker_scripts/pam-helper.sh 
COPY scripts/rsession.sh /rocker_scripts/rsession.sh 
COPY scripts/setup_R.sh /rocker_scripts/setup_R.sh 


RUN bash /rocker_scripts/install_rstudio.sh
RUN bash /rocker_scripts/install_pandoc.sh
RUN bash /rocker_scripts/install_quarto.sh

# Custom for RIS
RUN echo "server-user=${DEFAULT_USER}" >> /etc/rstudio/rserver.conf

# Extra packages for R
RUN sudo apt-get update && \
    sudo apt-get -y install libdeflate-dev libbz2-dev liblzma-dev libpcre2-dev libicu-dev libz-dev libbz2-dev liblzma-dev libxml2 libxml2-dev libxt6 libglpk-dev libmysqlclient21 gzip less wget gawk grep libnss-sss && \
    sudo apt-get clean

# Install conda (mamba version) with Jupyter
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/23.1.0-1/Mambaforge-23.1.0-1-Linux-x86_64.sh" && \
    bash Mambaforge-23.1.0-1-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Miniconda3-py310_23.3.1-0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/mamba install -y -c conda-forge jupyterlab nb_conda_kernels jupyter_contrib_nbextensions jupyterlab-spellchecker Jinja2 && \
    $CONDA_DIR/bin/pip3 install --root-user-action=ignore jupyterlab-quarto

# Create "basic-tools" and "gosling" conda environment and install some basic tools
COPY env/*.yml /tmp/
RUN $CONDA_DIR/bin/mamba env create -f /tmp/basic-tools.yml && \
    $CONDA_DIR/bin/mamba env create -f /tmp/hictools.yml && \
    $CONDA_DIR/bin/conda clean -y --all

# soft link conda bin
COPY scripts/conda_soft_links.sh /tmp/
RUN bash /tmp/conda_soft_links.sh $CONDA_DIR/bin

RUN Rscript -e 'install.packages(c("tidyverse","ggplot2","reshape2","cowplot","reticulate"))'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(version = "3.17")'
RUN Rscript -e 'BiocManager::install(c("bsseq","DESeq2","rtracklayer"))'

EXPOSE 8787

CMD ["/init"]
