FROM ubuntu:xenial
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL description="Heavy container for Chromoseq"

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
    python-dev \
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
    default-jdk \
    default-jre \
    g++ \
    git \
    less \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    pkg-config \
    python \
    unzip \
    zip \
    libbz2-dev \
    ca-certificates \
    file \
    g++ \
    libbz2-1.0 \
    libcurl3 \
    libicu55 \
    libopenblas-dev \
    liblzma5 \
    locales \
    zlib1g \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
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
		       

##############
#HTSlib 1.9#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.9.tar.bz2 && \
    cd /tmp/htslib-1.9 && \
    ./configure --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

################
#Samtools 1.9#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar --bzip2 -xf samtools-1.9.tar.bz2 && \
    cd /tmp/samtools-1.9 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.9


##############
## bedtools ##

WORKDIR /usr/local
RUN git clone https://github.com/arq5x/bedtools2.git && \
    cd /usr/local/bedtools2 && \
    git checkout v2.27.0 && \
    make && \
    ln -s /usr/local/bedtools2/bin/* /usr/local/bin/

############################
# R, bioconductor packages #
ARG R_VERSION
ENV R_VERSION=${R_VERSION:-3.6.0}
RUN cd /tmp/ && \
    ## Download source code
    curl -O https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz && \
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
   
##########################################
# Install conda and all python stuff
##########################################

# Configure environment
ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH

RUN cd /tmp && \
    mkdir -p $CONDA_DIR && \
    curl -s https://repo.continuum.io/miniconda/Miniconda3-4.3.21-Linux-x86_64.sh -o miniconda.sh && \
    /bin/bash miniconda.sh -f -b -p $CONDA_DIR && \
    rm miniconda.sh && \
    $CONDA_DIR/bin/conda config --system --add channels conda-forge && \
    $CONDA_DIR/bin/conda config --system --set auto_update_conda false && \
    conda clean -tipsy

RUN conda config --add channels bioconda && \
    conda install -c conda-forge petl && \
    conda install -c anaconda biopython scipy cython cyvcf2 && \
    conda install -y -c bioconda mosdepth

RUN cd /tmp && git clone https://github.com/pysam-developers/pysam.git && \
    cd pysam && \
    export HTSLIB_LIBRARY_DIR=$HTSLIB_INSTALL_DIR/lib && \
    export HTSLIB_INCLUDE_DIR=$HTSLIB_INSTALL_DIR/include && \
    python setup.py install

# Install Python 2 
RUN conda create --quiet --yes -p $CONDA_DIR/envs/python2 python=2.7 'pip' && \
    conda clean -tipsy && \
    /bin/bash -c "source activate python2 && \
    conda install -c anaconda svtools && \
    source deactivate"

#
#  install manta
#
ENV manta_version 1.5.0
WORKDIR /opt/
RUN wget https://github.com/Illumina/manta/releases/download/v${manta_version}/manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    tar -jxvf manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    mv manta-${manta_version}.centos6_x86_64 /usr/local/src/manta

#
# install varscan
# 

ENV VARSCAN_INSTALL_DIR=/opt/varscan

WORKDIR $VARSCAN_INSTALL_DIR
RUN wget https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar && \
  ln -s VarScan.v2.4.2.jar VarScan.jar

#
# pindel
#

WORKDIR /opt
RUN wget https://github.com/genome/pindel/archive/v0.2.5b8.tar.gz && \
  tar -xzf v0.2.5b8.tar.gz

WORKDIR /opt/pindel-0.2.5b8
RUN ./INSTALL $HTSLIB_INSTALL_DIR

WORKDIR /
RUN ln -s /opt/pindel-0.2.5b8/pindel /usr/local/bin/pindel && \
    ln -s /opt/pindel-0.2.5b8/pindel2vcf /usr/local/bin/pindel2vcf

##################
# ucsc utilities #

RUN mkdir -p /tmp/ucsc && \
    cd /tmp/ucsc && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod ugo+x * && \
    mv * /usr/local/bin/

############
# ichorCNA #
############

WORKDIR /usr/local/bin
RUN git clone https://github.com/broadinstitute/ichorCNA.git
RUN Rscript -e "install.packages(c('plyr', 'optparse','BiocManager')); BiocManager::install(c('HMMcopy','GenomeInfoDb','GenomicRanges'))"
RUN R CMD INSTALL ichorCNA


#
# Minimap
#

RUN cd /opt/ && git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make && \
    cp minimap2 /usr/local/bin/ && \
    rm -Rf minimap2


########
#VEP 90#
########

RUN cpan install DBI && cpan install Module::Build.pm && cpan install JSON 

RUN mkdir /opt/vep/
WORKDIR /opt/vep

RUN git clone https://github.com/Ensembl/ensembl-vep.git
WORKDIR /opt/vep/ensembl-vep
RUN git checkout postreleasefix/90

RUN perl INSTALL.pl --NO_UPDATE

WORKDIR /
RUN ln -s /opt/vep/ensembl-vep/vep /usr/bin/variant_effect_predictor.pl


#install docker, instructions from https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-using-the-repository
RUN apt-get update && apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common

RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -

RUN add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

RUN apt-get update

RUN apt-get install -y docker-ce

WORKDIR /opt/
RUN wget https://github.com/broadinstitute/cromwell/releases/download/36/cromwell-36.jar

#
# Cleanup
#

## Clean up
RUN cd / && \
   rm -rf /tmp/* && \
   apt-get autoremove -y && \
   apt-get autoclean -y && \
   rm -rf /var/lib/apt/lists/* && \
   apt-get clean && \
   rm -f /opt/*.bz2 /opt/*.gz

RUN mkdir /opt/files/

COPY addReadCountsToVcfCRAM.py /usr/local/bin/addReadCountsToVcfCRAM.py
COPY duphold_static /usr/local/bin/duphold_static
COPY FilterManta.pl /usr/local/bin/FilterManta.pl
COPY ichorToVCF.pl /usr/local/bin/ichorToVCF.pl
COPY make_report.py /usr/local/bin/make_report.py
COPY configManta.hg38.py.ini /opt/files/configManta.hg38.py.ini
COPY nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds /opt/files/nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds
COPY basespace_cromwell.config /opt/files/basespace_cromwell.config

COPY all_sequences.fa.bed.gz /opt/files/all_sequences.fa.bed.gz
COPY all_sequences.fa.bed.gz.tbi /opt/files/all_sequences.fa.bed.gz.tbi
COPY all_sequences.fa.fai /opt/files/all_sequences.fa.fai

COPY chromoseq_genes.bed /opt/files/chromoseq_genes.bed
COPY hg38.cytoBandIdeo.bed.gz /opt/files/hg38.cytoBandIdeo.bed.gz
COPY hg38.cytoBandIdeo.bed.gz.tbi /opt/files/hg38.cytoBandIdeo.bed.gz.tbi
COPY chromoseq_sv_filter.bedpe.gz /opt/files/chromoseq_sv_filter.bedpe.gz
COPY chromoseq_translocations.bedpe /opt/files/chromoseq_translocations.bedpe

COPY driver.py /opt/files/driver.py

RUN chmod a+wrx /opt/files/*
RUN chmod a+wrx /usr/local/bin/*

RUN rm -f /usr/local/lib/R/site-library/ichorCNA/extdata/*.rds /usr/local/lib/R/site-library/ichorCNA/extdata/MBC* /usr/local/lib/R/site-library/ichorCNA/extdata/*10kb* /usr/local/lib/R/site-library/ichorCNA/extdata/*50kb* /usr/local/lib/R/site-library/ichorCNA/extdata/*1000kb*
