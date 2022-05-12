# docker-regtools

FROM dhspence/docker-baseimage
MAINTAINER "Dave Spencer" <dhspence@gmail.com>
RUN apt-get update -y && \
    apt-get install libncurses5-dev cmake -y

RUN cd /opt/ && \
    git clone https://github.com/griffithlab/regtools && \
    cd regtools/ && \
    mkdir build && \
    cd build/ && \
    cmake .. && \
    make && \
    cp regtools /usr/local/bin/






