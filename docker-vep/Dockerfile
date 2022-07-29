FROM dhspence/docker-genomic-analysis:latest
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL description="Container with VEP"


########
#VEP 93#
########

RUN cpan install DBI && cpan install Module::Build.pm

RUN mkdir /opt/vep/
WORKDIR /opt/vep

RUN git clone https://github.com/Ensembl/ensembl-vep.git
WORKDIR /opt/vep/ensembl-vep
RUN git checkout postreleasefix/93

RUN perl INSTALL.pl --NO_UPDATE

RUN mkdir -p /opt/lib/perl/VEP/Plugins && chmod a+wrx /opt/lib/perl/VEP/Plugins && \
    git clone https://github.com/Ensembl/VEP_plugins.git && \
    mv VEP_plugins/* /opt/lib/perl/VEP/Plugins/

COPY Wildtype.pm /opt/lib/perl/VEP/Plugins/Wildtype.pm

WORKDIR /
RUN ln -s /opt/vep/ensembl-vep/vep /usr/bin/variant_effect_predictor.pl

COPY add_annotations_to_table_helper.py /usr/local/bin/add_annotations_to_table_helper.py

RUN chmod a+wrx /usr/local/bin/add_annotations_to_table_helper.py
