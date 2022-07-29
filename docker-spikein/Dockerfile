FROM dhspence/docker-baseimage
MAINTAINER "Dave Spencer" <dhspence@gmail.com>

RUN apt-get install unzip

WORKDIR /opt/

RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make && \
    cp minimap2 /usr/local/bin/    

# pear
ADD pear-0.9.6-bin-64 /usr/local/bin/
RUN mv /usr/local/bin/pear-0.9.6-bin-64 /usr/local/bin/pear && chmod a+x /usr/local/bin/pear

# minimap2
RUN wget https://github.com/ComputationalSystemsBiology/EoulsanDockerFiles/raw/master/TrimAdapt/cutadapt-1.8.1.tar.gz && \
    tar -xzf cutadapt-1.8.1.tar.gz && \
    cp /opt/cutadapt-1.8.1/bin/cutadapt /usr/local/bin/

# trimgalore
RUN wget https://github.com/ComputationalSystemsBiology/EoulsanDockerFiles/raw/master/TrimAdapt/trim_galore_v0.4.1.zip && \
    unzip trim_galore_v0.4.1.zip -d . && \
    cp /opt/trim_galore_zip/trim_galore /usr/local/bin/

# perl summarization script
ADD ProcessSpikein.pl /usr/local/bin/
