FROM ontresearch/dorado@sha256:sha256:1a98a3d8a5341d98bd84dd41ce2bc2a5a2dae5cb31268d0b3aff9dfec52defd6
LABEL description="Image with dorado and all models"
USER root

RUN export CUDA_VISIBLE_DEVICES=0

# Install essential packages
RUN sudo mkdir -p /opt/dorado/models && \
    dorado download --directory /opt/dorado/models --model all

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin