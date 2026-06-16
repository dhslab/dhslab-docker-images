FROM dhspence/docker-baseimage:latest

LABEL version="1.0"
LABEL description="DAF-seq tools: minimap2 + python deps for the nf-dafseq dedup & phasing steps"

ENV DEBIAN_FRONTEND=noninteractive

# The base image (python 3.8) already provides samtools, pysam, pandas, numpy, scipy,
# scikit-learn, matplotlib and the UCSC bedGraphToBigWig used by phasing. Add what it lacks:

# minimap2 (alignment step), pinned to match the workflow
ENV MINIMAP2_VERSION=2.30
RUN cd /opt && \
    git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    git checkout v${MINIMAP2_VERSION} && \
    make && \
    cp minimap2 /usr/local/bin/ && \
    cd /opt && rm -rf minimap2/*.o

# seaborn is the only QC-plotting dep missing from the base; install it without pinning.
# (Do NOT pin matplotlib here: the base image is python 3.8, and matplotlib>=3.6 drops py3.8 —
# pip would fail. scikit-learn and matplotlib are already present in the base.)
RUN pip install --no-cache-dir seaborn

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Sanity check: fail the build if anything the pipeline needs is missing
RUN minimap2 --version && \
    samtools --version | head -1 && \
    bedGraphToBigWig 2>&1 | head -1; \
    python3 -c "import pysam, pandas, numpy, scipy, sklearn, matplotlib, seaborn; print('dafseq py deps OK')"
