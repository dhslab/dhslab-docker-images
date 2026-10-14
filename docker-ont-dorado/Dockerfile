FROM ontresearch/dorado@sha256:832c995858ac5418f4357d87c31faf8259e4dba45201c78216704baa05d651c6
LABEL description="Image with dorado and all models"
USER root

RUN export CUDA_VISIBLE_DEVICES=0

# Install essential packages
RUN sudo mkdir -p /opt/dorado/models && \
    dorado download --directory /opt/dorado/models --model all

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
