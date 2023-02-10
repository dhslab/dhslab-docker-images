FROM ubuntu:20.04

MAINTAINER dhspence

RUN apt-get update --fix-missing && \
  apt-get install -q -y wget curl bzip2 \
  libbz2-dev git build-essential zlib1g-dev locales vim fontconfig ttf-dejavu


# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# Install conda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /usr/local -u -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Install conda dependencies
ADD environment.yml /

RUN pwd
RUN conda config --set always_yes yes --set changeps1 no && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --get && \
    conda update -q conda && \
    conda info -a && \
    conda install -c conda-forge ncurses && \
    conda env update -q -n root --file environment.yml && \
    conda clean --tarballs --index-cache

WORKDIR /tmp
RUN conda install tabulate && git clone https://github.com/dovetail-genomics/Micro-C.git && \
cp Micro-C/get_qc.py /usr/local/bin && \
rm -r /tmp/Micro-C

RUN conda install -c bioconda preseq bedtools>=2.28.0 mosdepth