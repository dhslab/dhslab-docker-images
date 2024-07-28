FROM ghcr.io/dhslab/docker-python3:latest

LABEL Image with python, htslib, samtools and whatshap

# update as needed
RUN apt-get update -y && apt-get install -y --no-install-recommends wget

##############
#HTSlib 1.18#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.18.tar.bz2 && \
    cd /tmp/htslib-1.18 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix && \
    ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip && \
    rm -Rf /tmp/htslib-1.18

################
#Samtools 1.18#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar --bzip2 -xf samtools-1.18.tar.bz2 && \
    cd /tmp/samtools-1.18 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.18


RUN pip3 install whatshap
