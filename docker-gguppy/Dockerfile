FROM nvidia/cuda:11.6.0-base-ubuntu20.04

# File Author / Maintainer
MAINTAINER David Spencer <dspencerlab@gmail.com>

LABEL Image with Guppy ONT basecaller and tools

ARG PACKAGE_VERSION=6.4.2
ARG BUILD_PACKAGES="wget apt-transport-https tzdata"
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update && \
    apt install --yes $BUILD_PACKAGES libnvidia-compute-460-server && \ 
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${PACKAGE_VERSION}-1~bionic_amd64.deb && \
    apt install --yes /tmp/ont_guppy_${PACKAGE_VERSION}-1~bionic_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata
