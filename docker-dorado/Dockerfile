FROM ontresearch/dorado@sha256:e70f10b2c11e8b772edf00d07107e8dec8da281a49e90f8908aa68085959a34d
LABEL description="Image with dorado"
USER root
# Install essential packages
RUN sudo mkdir -p /opt/dorado/models && \
    dorado download --directory /opt/dorado/models --model all
    
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
