FROM ubuntu:20.04
MAINTAINER Dave Spencer <dhspence@gmail.com>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y r-base-core
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('PSCBS', 'cghFLasso')"

RUN apt-get install -y \
    liblzma-dev \
    python-biopython \
    python-dev \
    python-matplotlib \
    python-numpy \
    python-pip \
    python-reportlab \
    python-scipy \
    python-tk \
    zlib1g-dev \
    emacs \
    emacs-goodies-el
    
RUN pip install -U future futures pandas pyfaidx pysam
RUN pip install cnvkit==0.9.1

# Let matplotlib build its font cache
RUN cnvkit.py version

# ENTRYPOINT ["cnvkit.py"]
# CMD ["--help"]

CMD ["bash"]
