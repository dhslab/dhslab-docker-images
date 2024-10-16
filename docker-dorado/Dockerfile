FROM ontresearch/dorado@sha256:2c9f8b2b21966f16479e66569c506fc81eedc31db7719f88892544857103dcbb
LABEL description="Image with dorado and all models"
USER root
# Install essential packages
RUN sudo mkdir -p /opt/dorado/models && \
    dorado download --directory /opt/dorado/models --model all
    
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
