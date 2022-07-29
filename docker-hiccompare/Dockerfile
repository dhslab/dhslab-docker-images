FROM dhspence/docker-genomic-analysis:122121
MAINTAINER David Spencer <dhspence@gmail.com>

LABEL description="hiccompare"

## install r packages, bioconductor, etc ##
ADD rpackages.R /tmp/
RUN R -f /tmp/rpackages.R

WORKDIR /opt

RUN git clone https://github.com/ay-lab/selfish
