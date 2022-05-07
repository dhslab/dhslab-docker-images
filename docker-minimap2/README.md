# docker-minimap2

FROM dhspence/docker-baseimage
MAINTAINER "Dave Spencer" <dhspence@gmail.com>
RUN apt-get update -y && \
    apt-get install libncurses5-dev -y

RUN cd /opt/ && \
    git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make && \
    cp minimap2 /usr/local/bin/

ADD pear-0.9.6-bin-64 /usr/local/bin/

RUN cp /usr/local/bin/pear-0.9.6-bin-64 pear

