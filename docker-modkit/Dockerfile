FROM ubuntu:focal-20221130

LABEL description="Basic image with modkit"

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget

# build modkit
RUN mkdir /opt/modkit && \
    wget https://github.com/nanoporetech/modkit/releases/download/v0.3.1rc1/modkit_v0.3.1rc1_centos7_x86_64.tar.gz && \
    tar -xf modkit_v0.3.1rc1_centos7_x86_64.tar.gz -C /opt/modkit && \
    chmod +x /opt/modkit/dist/modkit && \
    ln -s /opt/modkit/dist/modkit /usr/local/bin/modkit

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
