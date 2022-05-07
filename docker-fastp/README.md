# docker-fastp

FROM dhspence/docker-baseimage:031822
MAINTAINER "Dave Spencer" <dhspence@gmail.com>

WORKDIR /opt/

RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv ./fastp /usr/local/bin/



