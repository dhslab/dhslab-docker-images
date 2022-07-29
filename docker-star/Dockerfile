FROM dhspence/docker-baseimage

# File Author / Maintainer
MAINTAINER David Spencer <dspencer@wustl.edu>

COPY GeneralAlignment.sh /usr/local/bin/

RUN chmod a+x /usr/local/bin/GeneralAlignment.sh

##############
# STAR
##############
ENV STAR_VERSION 2.6.1d

ENV PACKAGES gcc g++ make wget zlib1g-dev

RUN set -ex

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    cd /opt/ && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && \
    tar xzf ${STAR_VERSION}.tar.gz && \
    cd STAR-${STAR_VERSION}/source && \
    make STARstatic && \
    cp STAR /usr/local/bin/ && \
    cd /opt/ && \
    'rm' -rf STAR-${STAR_VERSION} && \
    apt-get --purge autoremove -y  ${PACKAGES}

#RUN apt-get update -y && apt-get install -y --no-install-recommends \
#    libdb4.8-dev \
#    libdb4.8++-dev
    
RUN cpan install DB_File && \
   cpan install URI::Escape && \
   cpan install Set::IntervalTree && \
   cpan install Carp::Assert && \
   cpan install JSON::XS && \
   cpan install PerlIO::gzip

WORKDIR /usr/local/
RUN git clone --recursive https://github.com/STAR-Fusion/STAR-Fusion.git

RUN ln -s /usr/local/STAR-Fusion/STAR-Fusion /usr/local/bin/
