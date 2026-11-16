FROM mdivr/r-ver_4.2.2:v2

RUN mkdir -p /opt/chicago/bin

COPY chicago/chicagoTools/* /opt/chicago/bin/

RUN chmod -R 777 /opt/chicago/bin && \
    ln -s /opt/chicago/bin/* /usr/local/bin/ && \
    rm /usr/local/bin/*.md && \
    ln -s /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r && \
    Rscript -e "install.packages('BiocManager')" && \
    Rscript -e "BiocManager::install(version = '3.16')" && \
    installBioc.r --update --deps TRUE --error GenomicInteractions AnnotationHub && \
    install2.r --deps TRUE --error argparser devtools && \
    Rscript -e "library(devtools); install_bitbucket('chicagoTeam/Chicago', subdir='Chicago')"

# Install Chicdiff dependencies
RUN Rscript -e "BiocManager::install(c('DESeq2', 'IHW'), ask=FALSE, update=FALSE)"

# Install Chicdiff and ChicdiffData
RUN Rscript -e "remotes::install_github('RegulatoryGenomicsGroup/chicdiff', subdir='Chicdiff', upgrade='never')" && \
    Rscript -e "remotes::install_github('RegulatoryGenomicsGroup/chicdiff', subdir='ChicdiffData', upgrade='never')"

# Verify packages
RUN Rscript -e "stopifnot(requireNamespace('Chicago'), requireNamespace('Chicdiff'), requireNamespace('ChicdiffData'))"

CMD ["/bin/bash"]
