FROM mdivr/r-ver_4.2.2:v2

# Headers and shared libraries required when R dependencies build from source.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        libcurl4-openssl-dev \
        libssl-dev \
        libpcre2-dev \
        libuv1 \
        libuv1-dev \
        pkg-config \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/chicago/bin

COPY chicago/chicagoTools/* /opt/chicago/bin/

ENV R_DEFAULT_INTERNET_TIMEOUT=600

RUN Rscript -e "stopifnot(getOption('timeout') == 600)"

RUN chmod -R 777 /opt/chicago/bin && \
    ln -s /opt/chicago/bin/* /usr/local/bin/ && \
    rm /usr/local/bin/*.md && \
    ln -s /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r && \
    Rscript -e "install.packages('BiocManager'); stopifnot(requireNamespace('BiocManager'))"

# Keep the Bioconductor release compatible with R 4.2 without updating base-image packages.
RUN Rscript -e "BiocManager::install(version='3.16', ask=FALSE, update=FALSE)"

RUN for attempt in 1 2 3; do \
        Rscript -e "install.packages(c('curl', 'openssl', 'httr'), repos='https://cloud.r-project.org'); stopifnot(requireNamespace('curl'), requireNamespace('openssl'), requireNamespace('httr'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

# Retry interrupted downloads. Installed dependencies are reused on the next attempt.
# Omitting --update avoids rebuilding unrelated packages such as littler.
RUN for attempt in 1 2 3; do \
        installBioc.r --deps NA --error --skipinstalled GenomicInteractions AnnotationHub && \
        Rscript -e "stopifnot(requireNamespace('GenomicInteractions'), requireNamespace('AnnotationHub'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

RUN for attempt in 1 2 3; do \
        install2.r --deps NA --error argparser remotes && \
        Rscript -e "stopifnot(requireNamespace('argparser'), requireNamespace('remotes'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

RUN for attempt in 1 2 3; do \
        Rscript -e "remotes::install_bitbucket('chicagoTeam/Chicago', subdir='Chicago', dependencies=NA, upgrade='never'); stopifnot(requireNamespace('Chicago'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

# Install Chicdiff dependencies without optional suggested packages.
RUN for attempt in 1 2 3; do \
        Rscript -e "BiocManager::install(c('DESeq2', 'IHW'), dependencies=NA, ask=FALSE, update=FALSE); stopifnot(requireNamespace('DESeq2'), requireNamespace('IHW'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

RUN for attempt in 1 2 3; do \
        Rscript -e "options(repos=BiocManager::repositories()); remotes::install_github('RegulatoryGenomicsGroup/chicdiff', subdir='Chicdiff', dependencies=NA, upgrade='never'); stopifnot(requireNamespace('Chicdiff'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

RUN for attempt in 1 2 3; do \
        Rscript -e "remotes::install_github('RegulatoryGenomicsGroup/chicdiff', subdir='ChicdiffData', dependencies=NA, upgrade='never'); stopifnot(requireNamespace('ChicdiffData'))" && break; \
        if [ "$attempt" -eq 3 ]; then exit 1; fi; \
        sleep 10; \
    done

# Verify the installed packages in a fresh R process.
RUN Rscript -e "stopifnot(requireNamespace('Chicago'), requireNamespace('Chicdiff'), requireNamespace('ChicdiffData'), requireNamespace('GenomicInteractions'), requireNamespace('AnnotationHub'))"

CMD ["/bin/bash"]
