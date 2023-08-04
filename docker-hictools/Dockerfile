FROM ubuntu:20.04

LABEL Image for micro-c/hi-c tools

ENV DEBIAN_FRONTEND=noninteractive

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
    tzdata

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

## Conda installation
# Configure environment
ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH

# Install pythong/conda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh && \
    bash Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -p $CONDA_DIR -u -b && \
    rm Miniconda3-py37_23.1.0-1-Linux-x86_64.sh && \
    conda clean -y --all && \
    conda install -c bioconda samtools ucsc-bedgraphtobigwig && \ 
    pip install pysam tabulate numpy scipy py2bit matplotlib pandas h5py pyBigWig pairtools cooler fanc && \
    pip install cooltools && \
    ln -s ${CONDA_DIR}/bin/* /usr/local/bin/

RUN git clone https://github.com/billgreenwald/pgltools.git && mv pgltools /usr/local/ && \
    ln -s /usr/local/pgltools/sh/pgltools /usr/local/bin/pgltools && \
    sed -i -e 's!$script_path!/usr/local/pgltools/sh!g' /usr/local/pgltools/sh/pgltools

RUN pip install hicexplorer
