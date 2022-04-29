FROM dhspence/docker-genomic-analysis:latest
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL description="Heavy container for Chromoseq"

#
#  install manta
#
ENV manta_version 1.5.0
WORKDIR /opt/
RUN wget https://github.com/Illumina/manta/releases/download/v${manta_version}/manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    tar -jxvf manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    mv manta-${manta_version}.centos6_x86_64 /usr/local/src/manta

#
# install hmmcopy and ichor
# 
RUN apt-get update && \
    apt-get install -y build-essential \
    	    	       cmake \
		       python-dev \
    		       python-pip \
                       git \
                       wget \
                       autoconf \
                       zlib1g-dev \
  		       fort77 \
		       liblzma-dev  \
		       libblas-dev \
		       gfortran \
		       gcc-multilib \
		       gobjc++ \
		       aptitude \
		       libreadline-dev \
		       python-dev \
		       libpcre3 \
		       libpcre3-dev
		       
RUN cd /opt/ && \
    git config --global http.sslVerify false && \
    git clone --recursive https://github.com/shahcompbio/hmmcopy_utils.git && \
    cd /opt/hmmcopy_utils && \
    cmake . && \
    make && \
    cp bin/* /usr/local/bin/

RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('HMMcopy'); biocLite('GenomeInfoDb'); install.packages(c('devtools','optparse'))"
RUN Rscript --default-packages=devtools -e "install_github('broadinstitute/ichorCNA')"
RUN cd /opt/ && wget https://github.com/broadinstitute/ichorCNA/archive/master.zip && \
       unzip master.zip && mv ichorCNA-master/scripts/*.R /usr/local/bin/ && rm -Rf master.zip ichorCNA-master

ENV VARSCAN_INSTALL_DIR=/opt/varscan

WORKDIR $VARSCAN_INSTALL_DIR
RUN wget https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar && \
  ln -s VarScan.v2.4.2.jar VarScan.jar

#
# pindel
#
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && \
  tar xvjf samtools-1.2.tar.bz2

WORKDIR /opt/samtools-1.2
RUN make

WORKDIR /opt
RUN wget https://github.com/genome/pindel/archive/v0.2.5b8.tar.gz && \
  tar -xzf v0.2.5b8.tar.gz

WORKDIR /opt/pindel-0.2.5b8
RUN ./INSTALL /opt/samtools-1.2/htslib-1.2.1

WORKDIR /
RUN ln -s /opt/pindel-0.2.5b8/pindel /usr/local/bin/pindel && \
    ln -s /opt/pindel-0.2.5b8/pindel2vcf /usr/local/bin/pindel2vcf

###############
# Octopus
###############

#RUN conda install -y -c conda-forge -c bioconda octopus

#
# GATK
#

#GATK 3.6#
ENV maven_package_name apache-maven-3.3.9
ENV gatk_dir_name gatk-protected
ENV gatk_version 3.6
RUN cd /tmp/ && wget -q http://mirror.nohup.it/apache/maven/maven-3/3.3.9/binaries/apache-maven-3.3.9-bin.zip

# LSF: Comment out the oracle.jrockit.jfr.StringConstantPool.
RUN cd /tmp/ \
    && git clone --recursive https://github.com/broadgsa/gatk-protected.git \
    && cd /tmp/gatk-protected && git checkout tags/${gatk_version} \
    && sed -i 's/^import oracle.jrockit.jfr.StringConstantPool;/\/\/import oracle.jrockit.jfr.StringConstantPool;/' ./public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/varianteval/VariantEval.java \
    && mv /tmp/gatk-protected /opt/${gatk_dir_name}-${gatk_version}
RUN cd /opt/ && unzip /tmp/${maven_package_name}-bin.zip \
    && rm -rf /tmp/${maven_package_name}-bin.zip LICENSE NOTICE README.txt \
    && cd /opt/ \
    && cd /opt/${gatk_dir_name}-${gatk_version} && /opt/${maven_package_name}/bin/mvn verify -P\!queue \
    && mv /opt/${gatk_dir_name}-${gatk_version}/protected/gatk-package-distribution/target/gatk-package-distribution-${gatk_version}.jar /opt/GenomeAnalysisTK.jar \
    && rm -rf /opt/${gatk_dir_name}-${gatk_version} /opt/${maven_package_name}

#
# pysam and scipy
#

RUN conda config --add channels r && conda config --add channels bioconda && conda install -c conda-forge petl && \
    conda install -c anaconda biopython && conda install -c anaconda scipy

RUN conda install pysam

#
# blat
#

WORKDIR /usr/local/bin/
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod a+x blat

########
#VEP 90#
########

RUN cpan install DBI && cpan install Module::Build.pm

RUN mkdir /opt/vep/
WORKDIR /opt/vep

RUN git clone https://github.com/Ensembl/ensembl-vep.git
WORKDIR /opt/vep/ensembl-vep
RUN git checkout postreleasefix/90

RUN perl INSTALL.pl --NO_UPDATE

WORKDIR /
RUN ln -s /opt/vep/ensembl-vep/vep /usr/bin/variant_effect_predictor.pl


RUN conda install -y -c bioconda cyvcf2 htslib samtools deeptools mosdepth

RUN export PATH=$PATH:/opt/conda/bin/ && \
    /bin/bash -c "source activate python2 && conda install -c bioconda svtools && source deactivate"

#
# Cleanup
#

## Clean up
RUN cd / && \
   rm -rf /tmp/* && \
   apt-get autoremove -y && \
   apt-get autoclean -y && \
   rm -rf /var/lib/apt/lists/* && \
   apt-get clean && \
   rm -f /opt/*.bz2 /opt/*.gz
   

RUN mkdir -p /opt/lib/perl/VEP/Plugins && chmod a+wrx /opt/lib/perl/VEP/Plugins 
COPY Downstream.pm /opt/lib/perl/VEP/Plugins/Downstream.pm
COPY Wildtype.pm /opt/lib/perl/VEP/Plugins/Wildtype.pm

RUN mkdir /opt/files/

COPY add_annotations_to_table_helper.py /usr/local/bin/add_annotations_to_table_helper.py
COPY docm_and_coding_indel_selection.pl /usr/local/bin/docm_and_coding_indel_selection.pl
COPY runIchorCNA.R /usr/local/bin/runIchorCNA.R
COPY addReadCountsToVcfCRAM.py /usr/local/bin/addReadCountsToVcfCRAM.py
COPY configManta.hg38.py.ini /opt/files/configManta.hg38.py.ini
COPY ChromoSeq.hg38.bed /opt/files/ChromoSeq.hg38.bed
COPY GeneRegions.bed /opt/files/GeneRegions.bed
COPY ChromoSeq.translocations.fixed.v3.sorted.hg38.bedpe /opt/files/ChromoSeq.translocations.fixed.v3.sorted.hg38.bedpe
COPY ChromoSeqReporter.hg38.pl /usr/local/bin/ChromoSeqReporter.hg38.pl
COPY BlatContigs.pl /usr/local/bin/BlatContigs.pl
COPY pslScore.pl /usr/local/bin/pslScore.pl
COPY hg38.blacklist.merged.bed /opt/files/hg38.blacklist.merged.bed
COPY B38.callset.public.bedpe.gz /opt/files/B38.callset.public.bedpe.gz
COPY all.stranded.filtered.merged.bedpe.gz /opt/files/all.stranded.filtered.merged.bedpe.gz
COPY all.stranded.filtered.merged.bedpe.gz.tbi /opt/files/all.stranded.filtered.merged.bedpe.gz.tbi
COPY GeneCoverageRegions.bed /opt/files/GeneCoverageRegions.bed
COPY ChromoSeq.translocations.qc.bed /opt/files/ChromoSeq.translocations.qc.bed 
COPY nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.rds /opt/files/nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.rds

RUN chmod a+wrx /opt/files/*
RUN chmod a+wrx /usr/local/bin/*

#WORKDIR /opt/
#RUN conda install --yes 'pip' && \
#    conda clean -tipsy && \
#    pip install cmake==3.13
    
#RUN git clone -b master https://github.com/luntergroup/octopus.git && \
#    cd octopus && \
#    /opt/conda/bin/python ./scripts/install.py 
