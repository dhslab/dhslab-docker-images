# docker-fithic

FROM ubuntu:hirsute
MAINTAINER David Spencer <dhspence@gmail.com>

LABEL \
    description="Fithic"

ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update -y && apt-get install -y \
    ant \
    apt-utils \
    git \
    libncurses5-dev \
    ncurses-dev \
    unzip \
    wget \
    zlib1g-dev

RUN apt-get update -y && apt-get install -y python3.6 python3-pip python3-dev build-essential nodejs libbz2-dev liblzma-dev
RUN pip3 install --upgrade pip
RUN pip3 install biopython
RUN pip3 install pysam
RUN pip3 install scipy
RUN pip3 install statsmodels
RUN pip3 install fithic
RUN ln -s $(which python3) /usr/local/bin/python

COPY HiCKRy.py /usr/local/bin/
COPY makefiles.pl /usr/local/bin/
COPY juicer_tools.jar /usr/local/bin/
