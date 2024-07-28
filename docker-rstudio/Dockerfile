FROM rocker/rstudio:4.3.3

ARG DEFAULT_USER
COPY default_user.sh /rocker_scripts/default_user.sh
# Custom for RIS
RUN bash /rocker_scripts/default_user.sh "${DEFAULT_USER}" && \
    echo "server-user=${DEFAULT_USER}" >> /etc/rstudio/rserver.conf

# Extra packages for R
RUN sudo apt-get update && \
    sudo apt-get -y install libz-dev libbz2-dev liblzma-dev libxml2 libxml2-dev libxt6 libglpk-dev libmysqlclient21 gzip less wget gawk grep libnss-sss && \
    sudo apt-get clean

RUN Rscript -e 'install.packages(c("tidyverse","ggplot2","reshape2","cowplot","reticulate","dplyr","tidyr"))'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(version = "3.18")'
RUN Rscript -e 'BiocManager::install(c("bsseq","DESeq2","rtracklayer"))'

EXPOSE 8787

CMD ["/init"]
