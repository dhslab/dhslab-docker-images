FROM dhspence/docker-baseimage:latest
MAINTAINER "Dave Spencer" <dhspence@gmail.com>

RUN apt-get upgrade -y

RUN apt-get install --no-install-recommends -y git python python-dev python-pip bedtools && \
    pip install --upgrade pip && \
    pip install setuptools numpy

RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv bwa /usr/local/bin/
        
WORKDIR /opt/

RUN git clone --recursive https://github.com/aryeelab/guideseq.git && \
    cd guideseq && \
    pip install -r requirements.txt


