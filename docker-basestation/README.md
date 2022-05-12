# docker-basestation

FROM dhspence/docker-genomic-analysis:latest
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL \
  description="Updated supercontainer"

RUN cpan install Statistics::Basic && cpan install JSON && cpan install YAML::Tiny

RUN wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O /usr/local/bin/bs && chmod a+x /usr/local/bin/bs

ENV BASESPACE_ACCESS_TOKEN=d546e5027af4414b8a0da61a73a4aff4
ENV BASESPACE_API_SERVER=https://api.basespace.illumina.com


