FROM ontresearch/dorado@sha256:d5a340e6d5528ef61eee9cfcacbce01c6b7eefa1456c526f8630886e5087a7c5
LABEL description="Image with dorado and all models"
USER root

RUN export CUDA_VISIBLE_DEVICES=0

# Install essential packages
RUN sudo mkdir -p /opt/dorado/models && \
    dorado download --directory /opt/dorado/models --model all

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin