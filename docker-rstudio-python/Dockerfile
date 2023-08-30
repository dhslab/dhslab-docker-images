FROM dhspence/docker-rstudio:v040523

LABEL Rstudio image with python

RUN apt-get update -y && apt-get install -y --no-install-recommends ca-certificates \
    build-essential \
    zlib1g-dev \
    libcurl4 \
    liblz4-tool \
    make \
    less \
    curl \
    libcurl4-gnutls-dev \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bedtools \
    wget \
    libhdf5-dev \
    autotools-dev \
    automake \
    unzip \
    tzdata \ 
    libstdc++6

RUN Rscript -e 'install.packages("reticulate")'

RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(version = "3.16")'

RUN Rscript -e 'BiocManager::install(c("bsseq","DESeq2","rtracklayer"))'

RUN python -m ensurepip --upgrade
RUN pip install py2bit pyBigWig cooltools cooler PyGLtools fanc matplotlib cython pandas tabulate numpy scipy h5py

RUN export HTSLIB_CONFIGURE_OPTIONS=--enable-plugins && pip install pysam
