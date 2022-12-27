FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Image for nextflow"

# Install essential packages and Java
RUN apt-get update -y && \
    apt install wget -y && \
    apt install openjdk-17-jre -y && \
    apt-get clean

# Install Nextflow
RUN wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v22.10.4/nextflow-22.10.4-all | bash && \
    chmod 777 nextflow && \
    mv nextflow /usr/local/bin/nextflow

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin