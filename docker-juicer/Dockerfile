FROM mdivr/cuda-10.0:v1

RUN mkdir /opt/juicer/ && \
    wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar && \
    mv juicer_tools_1.22.01.jar /opt/juicer/

# install bwa
RUN git clone https://github.com/lh3/bwa.git && \
    mv bwa /opt/ && \
    cd /opt/bwa && \
    git checkout tags/v0.7.17 && \
    make && \
    chmod +x /opt/bwa/bwa && \
    ln -s /opt/bwa/bwa /usr/local/bin/

ENV PATH=/usr/local/cuda-10.0/bin:/usr/local/cuda-10.0/:$PATH

# docker build --no-cache -t mdivr/juicer_tools:v2 .
