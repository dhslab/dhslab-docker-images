FROM ursamajorlab/jammy-python:3.11
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Base Spencerlab image with basic tools for bioinformatic analyses

ENV DEBIAN_FRONTEND=noninteractive

#some basic tools
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    g++ \
    git \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    nodejs \
    pkg-config \
    rsync \
    unzip \
    wget \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bc \
    curl \
    tzdata \
    emacs \
    gawk \
    openssh-client \
    grep \
    evince \
    libnss-sss \
    bedtools \
    default-jdk \
    default-jre \
    git-lfs \
    smbclient \
    && apt-get clean all
    
RUN pip install numpy && \
    pip install cruzdb && \
    pip install cython && \
    pip install pyfaidx && \
    pip install cyvcf2 && \
    pip install pandas && \
    pip install scipy && \
    pip install scikit-learn && \
    pip install pysam && \
    pip install statsmodels && \
    pip install pyranges && \
    pip install pyyaml && \
    pip install gtfparse && \
    pip install globus-cli && \
    pip install nextflow && \
    pip install nf-core && \
    pip install boto3 && \
    pip install gget && \
    pip install openpyxl && \
    pip install SQLAlchemy && \
    pip install pymupdf

RUN wget https://github.com/fls-bioinformatics-core/GFFUtils/archive/refs/tags/v0.12.0.tar.gz && \
    tar xzf GFFUtils-0.10.3.tar.gz && \
    pip install -r ./GFFUtils-0.10.3/requirements.txt && \
    pip install ./GFFUtils-0.10.3/
    
#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

##############
#HTSlib 1.18#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.18.tar.bz2 && \
    cd /tmp/htslib-1.18 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix && \
    ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip && \
    rm -Rf /tmp/htslib-1.18

################
#Samtools 1.18#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar --bzip2 -xf samtools-1.18.tar.bz2 && \
    cd /tmp/samtools-1.18 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.18

################
#bcftools 1.18#
################
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 && \
    tar --bzip2 -xf bcftools-1.18.tar.bz2 && \
    cd /tmp/bcftools-1.18 && \
    make prefix=$BCFTOOLS_INSTALL_DIR && \
    make prefix=$BCFTOOLS_INSTALL_DIR install && \
    ln -s /opt/bcftools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/bcftools-1.18


##################
# ucsc utilities #
RUN mkdir -p /tmp/ucsc && \
    cd /tmp/ucsc && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    chmod ugo+x * && \
    mv * /usr/local/bin/

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && ./aws/install

RUN apt-get install -y --no-install-recommends groff

RUN git clone https://github.com/lh3/seqtk.git && cd seqtk && make && mv seqtk /usr/local/bin/

