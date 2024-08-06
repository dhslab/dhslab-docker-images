FROM ubuntu:focal-20221130

LABEL description="R 4.4.0 with Bioconductor 3.19 and core spencerlab packages"

RUN apt-get update && apt-get install -y --no-install-recommends locales && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    LC_ALL=en_US.UTF-8 && \
    LANG=en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    TERM=xterm

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    autoconf \
    zlib1g-dev \
    fort77 \
    liblzma-dev  \
    libblas-dev \
    gfortran \
    gcc-multilib \
    aptitude \
    libreadline-dev \
    libpcre3 \
    libpcre3-dev \
    bzip2 \
    curl \
    g++ \
    git \
    less \
    libcurl4 \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    pkg-config \
    unzip \
    zip \
    libbz2-dev \
    ca-certificates \
    file \
    g++ \
    libbz2-1.0 \
    libicu55 \
    libopenblas-dev \
    liblzma5 \
    locales \
    zlib1g \
    libbz2-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libicu-dev \
    libreadline-dev \
    liblzma-dev \
    libxt-dev \
    perl \
    xvfb \
    zlib1g-dev \
    libnss-sss

RUN apt-get update && \
    apt-get install -y 		       

############################
# R, bioconductor packages #
ARG R_VERSION
ENV R_VERSION=${R_VERSION:-4.4.0}
RUN cd /tmp/ && \
    ## Download source code
    curl -O https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz && \
    ## Extract source code
    tar -xf R-${R_VERSION}.tar.gz && \
    cd R-${R_VERSION} && \
    ## Set compiler flags
    R_PAPERSIZE=letter && \
    R_BATCHSAVE="--no-save --no-restore" && \
    R_BROWSER=xdg-open && \
    PAGER=/usr/bin/pager && \
    PERL=/usr/bin/perl && \
    R_UNZIPCMD=/usr/bin/unzip && \
    R_ZIPCMD=/usr/bin/zip && \
    R_PRINTCMD=/usr/bin/lpr && \
    LIBnn=lib && \
    AWK=/usr/bin/awk && \
    CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" && \
    CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" && \
    ## Configure options
    ./configure --enable-R-shlib \
               --enable-memory-profiling \
               --with-readline \
               --with-blas="-lopenblas" \
               --disable-nls \
               --without-recommended-packages && \
    ## Build and install
    make && \
    make install && \
    ## Add a default CRAN mirror
    echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site && \
    ## Add a library directory (for user-installed packages)
    mkdir -p /usr/local/lib/R/site-library && \
    chown root:staff /usr/local/lib/R/site-library && \
    chmod g+wx /usr/local/lib/R/site-library && \
    ## Fix library path
    echo "R_LIBS_USER='/usr/local/lib/R/site-library'" >> /usr/local/lib/R/etc/Renviron && \
    echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron

RUN Rscript -e 'install.packages(c("tidyverse","ggplot2","reshape2","cowplot","dplyr","tidyr","vcfR","reticulate"))'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(version = "3.19")'
RUN Rscript -e 'BiocManager::install(c("bsseq","DESeq2","rtracklayer","methylKit"))'
