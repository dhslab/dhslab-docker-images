FROM ghcr.io/dhslab/docker-rbase3.6.0

LABEL description="Container with ichorCNA"

############
# ichorCNA #
############

RUN Rscript -e ".libPaths(c('/usr/local/lib/R/site-library','/usr/local/lib/R/library')); install.packages(c('plyr', 'optparse','BiocManager', 'foreach', 'doMC'))"
RUN Rscript -e ".libPaths(c('/usr/local/lib/R/site-library')); BiocManager::install(c('HMMcopy','GenomeInfoDb','GenomicRanges'))"

RUN Rscript -e ".libPaths(c('/usr/local/lib/R/site-library')); BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"

COPY ichorCNA /opt/ichorCNA
RUN cp -r /opt/ichorCNA /usr/local/bin/ && \
    chown root:staff /usr/local/bin/ichorCNA && \
    chmod -R a+rwx /usr/local/bin/ichorCNA

WORKDIR /usr/local/bin
RUN R CMD INSTALL ichorCNA
