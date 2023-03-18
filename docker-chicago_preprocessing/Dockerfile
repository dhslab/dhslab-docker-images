FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Basic image with perl and bedtools"


# Install essential packages
RUN apt-get update -y && \
    apt-get install -y wget build-essential perl bedtools && \
    apt-get clean

# Add bam2chicago script
RUN mkdir -p /opt/chicago/bin
COPY bam2chicago.sh /opt/chicago/bin/
RUN chmod 777 /opt/chicago/bin/bam2chicago.sh && \
    ln -s /opt/chicago/bin/bam2chicago.sh /usr/local/bin/bam2chicago.sh

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin