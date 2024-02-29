FROM ubuntu:20.04

# NOTE(ariel) Install dependencies for HTSlib.
RUN apt update && apt install --yes --no-install-recommends \
  wget \
  autoconf \
  automake \
  make \
  clang \
  perl \
  git \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-gnutls-dev \
  libssl-dev \
  ca-certificates

# NOTE(ariel) Download and uncompress HTSlib tarball.
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
RUN tar xf htslib-1.17.tar.bz2

# NOTE(ariel) Build and install HTSlib.
RUN cd htslib-1.17 && ./configure && make && make install

# NOTE(ariel) Set environment variables to access HTSlib shared objects.
ENV LD_LIBRARY_PATH="/usr/local/lib"
ENV PATH="/usr/local/lib:$PATH"

# NOTE(ariel) Compile and install Haplotect.
#COPY . /haplotect
#RUN cd haplotect && bash ./compile.sh --release && cp ./haplotect /usr/local/bin

RUN git clone https://github.com/dhslab/chaplotect.git && \
    cd chaplotect && bash ./compile.sh --release && cp ./haplotect /usr/local/bin
