FROM mdivr/r-ver_4.2.2:v2

COPY chicago/chicagoTools/* /usr/local/bin/

RUN cp /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/ && \
    chmod +x /usr/local/bin/*.sh; chmod +x /usr/local/bin/*.py; chmod +x /usr/local/bin/*.r; chmod +x /usr/local/bin/*.R; chmod +x /usr/local/bin/*.pl && \
    Rscript -e "install.packages('BiocManager')" && \
    Rscript -e "BiocManager::install(version = '3.16')" && \
    installBioc.r --update --deps TRUE --error GenomicInteractions AnnotationHub && \
    install2.r --deps TRUE --error argparser devtools && \
    Rscript -e "library(devtools); install_bitbucket('chicagoTeam/Chicago', subdir='Chicago')" && \
    rm -rf /tmp/*


CMD ["bin/bash"]

# build  v1

# bsub -G compute-dspencer -q general-interactive -Is -a 'docker_build(mdivr/chicago:v1)' -- --tag mdivr/chicago:v1 --no-cache -f /storage1/fs1/dspencer/Active/spencerlab/mohamed/apps/docker/chicago/v1/Dockerfile /storage1/fs1/dspencer/Active/spencerlab/mohamed/apps/docker/chicago/v1