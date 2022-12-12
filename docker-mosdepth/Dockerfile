FROM ubuntu:20.04
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Basic image with conda and mosdepth

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends locales && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    LC_ALL=en_US.UTF-8 && \
    LANG=en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    TERM=xterm

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    wget \
    tzdata 

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
#LSF: Java bug that need to change the /etc/timezone.
#     The above /etc/localtime is not enough.
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

# needed for MGI data mounts
RUN apt-get update && apt-get install -y libnss-sss gawk openssh-client grep evince && apt-get clean all

ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH

RUN cd /tmp && mkdir -p $CONDA_DIR 
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh 
RUN /bin/bash Miniconda3-latest-Linux-x86_64.sh -f -b -p $CONDA_DIR

RUN rm -f Miniconda3-latest-Linux-x86_64.sh  && \
    $CONDA_DIR/bin/conda config --system --add channels conda-forge && \
    $CONDA_DIR/bin/conda config --system --set auto_update_conda false && \
    conda clean -tipsy

RUN conda update -n base -c defaults conda && \
    conda install -c bioconda mosdepth

