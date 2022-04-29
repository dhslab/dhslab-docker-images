FROM dhspence/docker-baseimage:031822

MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image for basic ad-hoc bioinformatic analyses

# install hmmcopy
RUN apt-get update && apt-get install cmake \
				      build-essential \
				      fort77 \
				      xorg-dev \
				      liblzma-dev  \
				      libblas-dev \
				      gfortran \
				      gcc-multilib \
				      gobjc++ \
				      aptitude \
				      libreadline-dev \
				      libpcre3 libpcre3-dev --no-install-recommends -y && \
    cd /opt/ && \
    git config --global http.sslVerify false && \
    git clone --recursive https://github.com/shahcompbio/hmmcopy_utils.git && \
    cd /opt/hmmcopy_utils && \
    cmake . && \
    make && \
    cp bin/* /usr/local/bin/

# install R
ARG R_VERSION
ARG BUILD_DATE
ENV BUILD_DATE 2017-06-20
ENV R_VERSION=${R_VERSION:-3.5.0}
RUN cd /tmp/ && \
    ## Download source code
    curl -O https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz && \
    ## Extract source code
    tar -xf R-${R_VERSION}.tar.gz && \
    cd R-${R_VERSION} && \
    ## Set compiler flags
    R_PAPERSIZE=letter && \
    R_BATCHSAVE="--no-save --no-restore" && \
    R_BROWSER=xdg-open && \
    PAGER=/usr/bin/pager && \
    PERL=/usr/bin/perl && \
    R_UNZIPCMD=/usr/bin/unzip && \
    R_ZIPCMD=/usr/bin/zip && \
    R_PRINTCMD=/usr/bin/lpr && \
    LIBnn=lib && \
    AWK=/usr/bin/awk && \
    CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" && \
    CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" && \
    ## Configure options
    ./configure --enable-R-shlib \
               --enable-memory-profiling \
               --with-readline \
               --with-blas="-lopenblas" \
               --disable-nls \
	       --with-x=no \
               --without-recommended-packages && \
    ## Build and install
    make && \
    make install && \
    ## Add a default CRAN mirror
    echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site && \
    ## Add a library directory (for user-installed packages)
    mkdir -p /usr/local/lib/R/site-library && \
    chown root:staff /usr/local/lib/R/site-library && \
    chmod g+wx /usr/local/lib/R/site-library && \
    ## Fix library path
    echo "R_LIBS_USER='/usr/local/lib/R/site-library'" >> /usr/local/lib/R/etc/Renviron && \
    echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron && \
    ## Use littler installation scripts
    Rscript -e "install.packages(c('littler', 'docopt'))" && \
    ln -s /usr/local/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r && \
    ln -s /usr/local/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r && \
    ln -s /usr/local/lib/R/site-library/littler/bin/r /usr/local/bin/r

RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('HMMcopy'); biocLite('GenomeInfoDb'); install.packages(c('devtools','optparse'))"
RUN Rscript --default-packages=devtools -e "install_github('broadinstitute/ichorCNA', '--no-docs')"

RUN cd /opt/ && wget https://github.com/broadinstitute/ichorCNA/archive/master.zip && \
       unzip master.zip && mv ichorCNA-master/scripts/*.R /usr/local/bin/ && rm -Rf master.zip ichorCNA-master
       
## Clean up
RUN cd / && \
    rm -rf /tmp/* && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean
