FROM dhspence/docker-genomic-analysis:latest
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL description="Heavy container with Manta"

ENV manta_version 1.5.0

WORKDIR /opt/

RUN wget https://github.com/Illumina/manta/releases/download/v${manta_version}/manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    tar -jxvf manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    mv manta-${manta_version}.centos6_x86_64 /usr/local/src/manta

RUN apt-get update && apt-get install cmake \
				      build-essential \
				      fort77 \
				      xorg-dev \
				      liblzma-dev  \
				      libblas-dev \
				      gfortran \
				      gcc-multilib \
				      gobjc++ \
				      aptitude \
				      libreadline-dev \
				      libpcre3 libpcre3-dev --no-install-recommends -y && \
    cd /opt/ && \
    git config --global http.sslVerify false && \
    git clone --recursive https://github.com/shahcompbio/hmmcopy_utils.git && \
    cd /opt/hmmcopy_utils && \
    cmake . && \
    make && \
    cp bin/* /usr/local/bin/


