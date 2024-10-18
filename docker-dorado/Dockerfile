FROM ontresearch/dorado@sha256:2c9f8b2b21966f16479e66569c506fc81eedc31db7719f88892544857103dcbb
LABEL description="Image with dorado and all models"
USER root
# Install essential packages
RUN sudo mkdir -p /opt/dorado/models && \
    dorado download --directory /opt/dorado/models --model all

##############
#HTSlib 1.21#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.21.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.21.tar.bz2 && \
    cd /tmp/htslib-1.21 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    cp $HTSLIB_INSTALL_DIR/bin/tabix /usr/local/bin/tabix && \
    cp $HTSLIB_INSTALL_DIR/bin/bgzip /usr/local/bin/bgzip && \
    rm -rf /tmp/htslib-1.21

################
#Samtools 1.21#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar --bzip2 -xf samtools-1.21.tar.bz2 && \
    cd /tmp/samtools-1.21 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.21

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
