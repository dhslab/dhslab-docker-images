FROM dhspence/docker-genomic-analysis
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image with sicer

RUN apt-get update && apt-get install -y --no-install-recommends locales && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    LC_ALL=en_US.UTF-8 && \
    LANG=en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    TERM=xterm

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    wget \
    tzdata 

# needed for MGI data mounts
RUN apt-get update && apt-get install -y libnss-sss gawk openssh-client grep evince && apt-get clean all

WORKDIR /opt/
RUN git clone https://github.com/dariober/SICERpy.git && ln -s /opt/SICERpy/SICERpy/SICER.py /usr/local/bin/

RUN conda install -c bioconda epic2
RUN export PATH=$PATH:/opt/conda/bin/
RUN /bin/bash -c "source activate python2 && conda install -c bioconda macs2 && source deactivate"

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
#LSF: Java bug that need to change the /etc/timezone.
#     The above /etc/localtime is not enough.
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

