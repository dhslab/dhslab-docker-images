# docker-whatshap

FROM dhspence/docker-baseimagewithconda
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image for with freebayes and whatshapp

#some basic tools
RUN apt-get update -y && apt-get install -y --no-install-recommends

RUN pip3 install whatshap

COPY freebayes /usr/local/bin/

RUN chmod a+x /usr/local/bin/freebayes




