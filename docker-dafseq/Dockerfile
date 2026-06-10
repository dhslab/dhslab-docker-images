FROM dhspence/docker-baseimage:latest

LABEL version="1.0"
LABEL description="DAF-seq tools: minimap2 + python deps for the nf-dafseq dedup & phasing steps"

ENV DEBIAN_FRONTEND=noninteractive

# The base image already provides samtools, pysam, pandas, numpy, scipy and the UCSC
# bedGraphToBigWig used by phasing. Add what it lacks for nf-dafseq:

# minimap2 (alignment step), pinned to match the workflow
ENV MINIMAP2_VERSION=2.30
RUN cd /opt && \
    git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    git checkout v${MINIMAP2_VERSION} && \
    make && \
    cp minimap2 /usr/local/bin/ && \
    cd /opt && rm -rf minimap2/*.o

# python packages used by the dedup (clustering) and QC plotting that are not in the base image
RUN pip install --no-cache-dir \
    scikit-learn \
    seaborn \
    matplotlib==3.10.3

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Sanity check: fail the build if anything the pipeline needs is missing
RUN minimap2 --version && \
    samtools --version | head -1 && \
    bedGraphToBigWig 2>&1 | head -1; \
    python3 -c "import pysam, pandas, numpy, scipy, sklearn, matplotlib, seaborn; print('dafseq py deps OK')"
