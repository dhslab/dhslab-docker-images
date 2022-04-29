FROM dhspence/docker-baseimage
MAINTAINER "Dave Spencer" <dhspence@gmail.com>

COPY subread-1.6.3-Linux-x86_64.tar.gz /opt/

WORKDIR /opt/

RUN tar -zxvf /opt/subread-1.6.3-Linux-x86_64.tar.gz && \
    cp /opt/subread-1.6.3-Linux-x86_64/bin/featureCounts /usr/local/bin/

