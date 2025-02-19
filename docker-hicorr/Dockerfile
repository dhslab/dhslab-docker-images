FROM ubuntu:20.04

LABEL maintainer="Andrew Bonney abonney@wustl.edu"
LABEL description="Docker image for HiCorr micro-C/Hi-C analysis"

ENV DEBIAN_FRONTEND=noninteractive

# Install basic system dependencies
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    g++ \
    git \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    make \
    perl \
    python3 \
    python3-dev \
    python3-pip \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bedtools \
    bowtie \
    samtools \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

# Install required Perl modules
RUN cd /opt/ \
    && curl https://raw.github.com/miyagawa/cpanminus/master/cpanm > cpanm \
    && chmod +x cpanm \
    && ln -s /opt/cpanm /usr/local/bin/ \
    && cpanm install Statistics::Basic \
    && cpanm install JSON \
    && cpanm install List::Util \
    && cpanm install YAML::Tiny

# Install Python packages
RUN pip3 install numpy pandas scipy

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Set up entrypoint script
COPY entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
