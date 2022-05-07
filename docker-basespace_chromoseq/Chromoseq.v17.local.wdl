workflow ChromoSeq {

  String Cram
  String CramIndex 
  String Name
  String Gender
  String MappingSummary
  String? CoverageSummary
  String TumorCounts
  String OutputDir
  
  String Translocations
  String GenesBed
  
  String Cytobands
  String SVDB

  String CustomAnnotationVcf #   = "/gscmnt/gc3042/cle_validation/myeloseq_haloplex/accessory_files/myeloseq_custom_annotations.annotated.011618.b37.vcf.gz"
  String CustomAnnotationIndex # = "/gscmnt/gc3042/cle_validation/myeloseq_haloplex/accessory_files/myeloseq_custom_annotations.annotated.011618.b37.vcf.gz.tbi"
  String CustomAnnotationParameters #= "MYELOSEQ,vcf,exact,0,TCGA_AC,MDS_AC,MYELOSEQBLACKLIST"

  String HotspotVCF
  String MantaConfig
  
  String Reference
  String ReferenceIndex
  String ReferenceBED
  String VEP

  String gcWig
  String mapWig
  String ponRds
  String centromeres
  String genomeStyle
  String genome

  String tmp
  
  Float minVarFreq = 0.02
  Int MinReads = 3
  Float varscanPvalindel = 0.1
  Float varscanPvalsnv = 0.01

  Int MinCNASize = 2000000
  Float MinCNAabund = 5.0
  Int LowCNASize = 50000000
  Float LowCNAabund = 10.0
  
  String JobGroup  

  String? chromoseq_docker

  call prepare_bed {
    input: Bedpe=Translocations,
    Bed=GenesBed,
    Reference=ReferenceBED,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
    
  call cov_qc as gene_qc {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Name=Name,
    Bed=GenesBed,
    refFasta=Reference,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }

  call cov_qc as sv_qc {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Name=Name,
    Bed=prepare_bed.svbed,
    refFasta=Reference,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
  
  call run_manta {
    input: Bam=Cram,
    BamIndex=CramIndex,
    Config=MantaConfig,
    Reference=Reference,
    ReferenceBED=ReferenceBED,
    Name=Name,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
    
  call run_ichor {
    input: Bam=Cram,
    BamIndex=CramIndex,
    refFasta=Reference,
    ReferenceBED=ReferenceBED,
    tumorCounts=TumorCounts,
    gender=Gender,
    gcWig=gcWig,
    mapWig=mapWig,
    ponRds=ponRds,
    centromeres=centromeres,
    Name=Name,
    genomeStyle = genomeStyle,
    genome = genome,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
  
  call run_varscan {
    input: Bam=Cram,
    BamIndex=CramIndex,
    CoverageBed=GenesBed,
    MinFreq=minVarFreq,
    pvalsnv=varscanPvalsnv,
    pvalindel=varscanPvalindel,
    refFasta=Reference,
    Name=Name,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
  
  call run_pindel_region as run_pindel_flt3itd {
    input: Bam=Cram,
    BamIndex=CramIndex,
    Reg='chr13:28033987-28034316',
    refFasta=Reference,
    Name=Name,
    genome=genome,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }

  call combine_variants {
    input: VCFs=[run_varscan.varscan_snv_file,
    run_varscan.varscan_indel_file,
    run_pindel_flt3itd.pindel_vcf_file,
    HotspotVCF],
    Bam=Cram,
    BamIndex=CramIndex,
    refFasta=Reference,
    Name=Name,
    MinReads=MinReads,
    MinVAF=minVarFreq,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
  
  call annotate_variants {
    input: Vcf=combine_variants.combined_vcf_file,
    refFasta=Reference,
    Vepcache=VEP,
    Cytobands=Cytobands,
    CustomAnnotationVcf=CustomAnnotationVcf,
    CustomAnnotationIndex=CustomAnnotationIndex,
    CustomAnnotationParameters=CustomAnnotationParameters,
    Name=Name,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
  
  call annotate_svs {
    input: Vcf=run_manta.vcf,
    CNV=run_ichor.seg,
    refFasta=Reference,
    refFastaIndex=ReferenceIndex,
    Vepcache=VEP,
    SVAnnot=SVDB,
    Translocations=Translocations,
    Cytobands=Cytobands,
    minCNAsize=MinCNASize,
    minCNAabund=MinCNAabund,
    lowCNAsize=LowCNASize,
    lowCNAabund=LowCNAabund,    
    Name=Name,
    gender=Gender,
    jobGroup=JobGroup,
    tmp=tmp,
    docker=chromoseq_docker
  }
  
  call make_report {
    input: SVVCF=annotate_svs.vcf,
    GeneVCF=annotate_variants.annotated_filtered_vcf,
    KnownGenes=prepare_bed.genes,
    GeneQC=gene_qc.qc_out,
    SVQC=sv_qc.qc_out,
    MappingSummary=MappingSummary,
    CoverageSummary=CoverageSummary,
    Name=Name,
    jobGroup=JobGroup,
    docker=chromoseq_docker,
    tmp=tmp
  }
  
#  call make_igv {
#    input: Name=Name
#  }
  
  call gather_files {
    input: OutputFiles=[annotate_svs.vcf,
    annotate_svs.vcf_index,
    run_ichor.params,
    run_ichor.seg,
    run_ichor.genomewide_pdf,
    run_ichor.allgenomewide_pdf,
    run_ichor.rdata,run_ichor.wig,
    run_ichor.correct_pdf,
    gene_qc.qc_out,
    gene_qc.region_dist,
    gene_qc.global_dist,
    sv_qc.qc_out,
    sv_qc.region_dist,
    annotate_variants.annotated_filtered_vcf,
    make_report.report],  #make_bw.bigwig_file,
#    make_igv.igv_xml],
    OutputDir=OutputDir,
    jobGroup=JobGroup,
    docker=chromoseq_docker
  }
  
}

task prepare_bed {
  String Bedpe
  String Bed
  String Reference
  String jobGroup
  String? tmp
  String? docker = "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
  
  command <<<
    awk -v OFS="\t" '{ split($7,a,"_"); print $1,$2,$3,a[1],".",$9; print $4,$5,$6,a[2],".",$10; }' ${Bedpe} | sort -u -k 1,1V -k 2,2n > sv.bed
    ((cat sv.bed | cut -f 4) && (cat ${Bed} | cut -f 6)) > genes.txt
    gunzip -c ${Reference} | cut -f 1 > chroms.txt
  >>>

  runtime {
    docker_image: docker
    cpu: "1"
    memory: "4 G"
    job_group: jobGroup
  }

  output {
    File svbed = "sv.bed"
    File genes = "genes.txt"
    Array[String] chroms = read_lines("chroms.txt")
  }
}

task cov_qc {
  String Cram
  String CramIndex
  String Bed
  String Name
  String refFasta
  String jobGroup
  String tmp
  String? docker = "dhspence/docker-mosdepth"
  
  command <<<
    set -eo pipefail && \
    /opt/conda/bin/mosdepth -n -f ${refFasta} -t 4 -i 2 -x -Q 20 -b ${Bed} --thresholds 10,20,30,40 "${Name}" ${Cram} && \
    /usr/local/bin/bedtools intersect -header -b "${Name}.regions.bed.gz" -a "${Name}.thresholds.bed.gz" -wo | \
    awk -v OFS="\t" '{ if (NR==1){ print $0,"%"$5,"%"$6,"%"$7,"%"$8,"MeanCov"; } else { print $1,$2,$3,$4,$5,$6,$7,$8,sprintf("%.2f\t%.2f\t%.2f\t%.2f",$5/$NF*100,$6/$NF*100,$7/$NF*100,$8/$NF*100),$(NF-1); } }' > "${Name}."$(basename ${Bed} .bed)".covqc.txt" && \
    mv "${Name}.mosdepth.region.dist.txt" "${Name}.mosdepth."$(basename ${Bed} .bed)".region.dist.txt"
  >>>
  
  runtime {
    docker_image: docker
    cpu: "4"
    memory: "32 G"
    job_group: jobGroup
  }
  
  output {
    File qc_out = glob("*.covqc.txt")[0]
    File global_dist = "${Name}.mosdepth.global.dist.txt"
    File region_dist = glob("*.region.dist.txt")[0]
  }

}

task run_manta {
  String Bam
  String BamIndex 
  String Config
  String Name
  String Reference
  String ReferenceBED
  String jobGroup
  String tmp
  String? docker = "mgibio/basespace_chromoseq:v12"
  
  command <<<
    set -eo pipefail && \
    /usr/local/src/manta/bin/configManta.py --config=${Config} --tumorBam=${Bam} --referenceFasta=${Reference} \
    --runDir=manta --callRegions=${ReferenceBED} --outputContig && \
    ./manta/runWorkflow.py -m local -q research-hpc -j 4 -g 32 && \
    zcat ./manta/results/variants/tumorSV.vcf.gz | /bin/sed 's/DUP:TANDEM/DUP/g' > fixed.vcf && \
    /usr/local/bin/duphold_static -v fixed.vcf -b ${Bam} -f ${Reference} -t 4 -o ${Name}.tumorSV.vcf && \
    bgzip ${Name}.tumorSV.vcf && /usr/bin/tabix ${Name}.tumorSV.vcf.gz
  >>>
  runtime {
    docker_image: docker
    cpu: "4"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File vcf = "${Name}.tumorSV.vcf.gz"
    File index = "${Name}.tumorSV.vcf.gz.tbi"
  }
}

task count_reads {
  String Bam
  String BamIndex
  String ReferenceBED
  String Chrom
  String jobGroup
  String refFasta
  String refIndex
  String tmp
  String? docker = "mgibio/basespace_chromoseq:v12"
  
  command {
    set -eo pipefail && \
    (/usr/local/bin/bedtools makewindows -b ${ReferenceBED} -w 500000 | \
    awk -v OFS="\t" -v C="${Chrom}" '$1==C && NF==3' > ${tmp}/${Chrom}.windows.bed) && \
    /usr/local/bin/samtools view -b -f 0x2 -F 0x400 -q 20 -T ${refFasta} ${Bam} ${Chrom} | \
    /usr/local/bin/intersectBed -sorted -nobuf -c -bed -b stdin -a ${tmp}/${Chrom}.windows.bed > ${Chrom}.counts.bed
  }

  runtime {
    docker_image: docker
    cpu: "1"
    memory: "8 G"
    job_group: jobGroup
  }
  output {
    File counts_bed = "${Chrom}.counts.bed"
  }
}


task run_ichor {
  String Bam
  String BamIndex
  String ReferenceBED
  String tumorCounts
  String refFasta
  String Name
  String gender
  String genome
  String genomeStyle
  String jobGroup
  String gcWig
  String mapWig
  String ponRds
  String centromeres
  
  Int? minCNAsize
  Float? lowAbundVal
  Int? lowAbundCnaSize
  
  String? tmp
  String? docker = "mgibio/basespace_chromoseq:v12"
  
  command <<<
    set -eo pipefail && \
    tail -n +6 ${tumorCounts} | sort -k 1V,1 -k 2n,2 | awk -v window=500000 'BEGIN { chr=""; } { if ($1!=chr){ printf("fixedStep chrom=%s start=1 step=%d span=%d\n",$1,window,window); chr=$1; } print $5; }' > "${Name}.tumor.wig" && \
    /usr/local/bin/Rscript /gscmnt/gc2555/spencer/dhs/git/ichorCNA/scripts/runIchorCNA.R --id ${Name} \
    --WIG "${Name}.tumor.wig" --ploidy "c(2)" --normal "c(0.1,0.5,.85)" --maxCN 3 \
    --gcWig ${gcWig} \
    --mapWig ${mapWig} \
    --centromere ${centromeres} \
    --normalPanel ${ponRds} \
    --genomeBuild ${genome} \
    --sex ${gender} \
    --includeHOMD False --chrs "c(1:22, \"X\", \"Y\")" --chrTrain "c(1:22)" --fracReadsInChrYForMale 0.0005 \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --txnE 0.999999 --txnStrength 1000000 --genomeStyle ${genomeStyle} --outDir ./ --libdir /gscmnt/gc2555/spencer/dhs/git/ichorCNA/ && \
    awk -v G=${gender} '$2!~/Y/ || G=="male"' "${Name}.seg.txt" > "${Name}.segs.txt" && \
    mv ${Name}/*.pdf .
  >>>
  
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "16 G"
    job_group: jobGroup
  }
  
  output {
    File params = "${Name}.params.txt"
    File seg = "${Name}.segs.txt"
    File genomewide_pdf = "${Name}_genomeWide.pdf"
    File allgenomewide_pdf = "${Name}_genomeWide_all_sols.pdf"
    File correct_pdf = "${Name}_genomeWideCorrection.pdf"
    File rdata = "${Name}.RData"
    File wig = "${Name}.tumor.wig"
  }
}

task run_varscan {
  String Bam
  String BamIndex
  Int? MinCov
  Float? MinFreq
  Int? MinReads
  Float? pvalindel
  Float? pvalsnv
  String CoverageBed
  String refFasta
  String Name
  String jobGroup
  String? tmp
  String? docker = "mgibio/basespace_chromoseq:v12"
  
  command <<<
    /usr/local/bin/samtools mpileup -f ${refFasta} -l ${CoverageBed} ${Bam} > ${tmp}/mpileup.out && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2snp ${tmp}/mpileup.out --min-coverage ${default=6 MinCov} --min-reads2 ${default=3 MinReads} \
    --min-var-freq ${default="0.02" MinFreq} --p-value ${default="0.01" pvalsnv} --output-vcf | bgzip -c > ${Name}.snv.vcf.gz && tabix ${Name}.snv.vcf.gz && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2indel ${tmp}/mpileup.out --min-coverage ${default=6 MinCov} --min-reads2 ${default=3 MinReads} \
    --min-var-freq ${default="0.02" MinFreq} --p-value ${default="0.1" pvalindel} --output-vcf | bgzip -c > ${Name}.indel.vcf.gz && tabix ${Name}.indel.vcf.gz
  >>>
  
  runtime {
    docker_image: docker
    cpu: "2"
    memory: "16 G"
    job_group: jobGroup
  }
  output {
    File varscan_snv_file = "${Name}.snv.vcf.gz"
    File varscan_indel_file = "${Name}.indel.vcf.gz"
  }
}

task run_pindel_region {
  String Bam
  String BamIndex
  String Reg
  Int? Isize
  Int? MinReads
  String refFasta
  String Name
  String jobGroup
  String? tmp
  String genome
  String? docker
  
  command <<<
    (set -eo pipefail && /usr/local/bin/samtools view -T ${refFasta} ${Bam} ${Reg} | /opt/pindel-0.2.5b8/sam2pindel - ${tmp}/in.pindel ${default=250 Isize} tumor 0 Illumina-PairEnd) && \
    /usr/local/bin/pindel -f ${refFasta} -p ${tmp}/in.pindel -c ${Reg} -o ${tmp}/out.pindel && \
    /usr/local/bin/pindel2vcf -P ${tmp}/out.pindel -G -r ${refFasta} -e ${default=3 MinReads} -R ${default="hg38" genome} -d ${default="hg38" genome} -v ${tmp}/pindel.vcf && \
    /bin/sed 's/END=[0-9]*\;//' ${tmp}/pindel.vcf | bgzip -c > ${Name}.pindel.vcf.gz && tabix ${Name}.pindel.vcf.gz
  >>>
  
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "16 G"
    job_group: jobGroup
  }
  output {
    File pindel_vcf_file = "${Name}.pindel.vcf.gz"
  }
}

task run_platypus {
  String Bam
  String BamIndex
  String CoverageBed
  String? DocmVcf
  Float? MinFreq
  String Name
  String refFasta
  String jobGroup
  String? tmp
  String? docker
  
  command <<<
    /usr/bin/awk '{ print $1":"$2+1"-"$3; }' ${CoverageBed} > "regions.txt" && \
    /opt/conda/bin/octopus -R ${refFasta} -I ${Bam} -t regions.txt -C cancer > "${Name}.vcf" && \
    /bin/sed 's/VCFv4.3/VCFv4.1/' "${Name}.vcf" > "${Name}.platypus.vcf"     
  >>>
  
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File platypus_vcf_file = "${Name}.platypus.vcf"
  }
}

task combine_variants {
  Array[String] VCFs
  String Bam
  String BamIndex
  String refFasta
  String Name
  Int MinReads
  Float MinVAF
  String jobGroup
  String? tmp
  String? docker = "mgibio/basespace_chromoseq:v12"

  command {
    /opt/conda/envs/python2/bin/bcftools merge --force-samples -O z ${sep=" " VCFs} | \
    /opt/conda/envs/python2/bin/bcftools norm -m- -f ${refFasta} -O z > ${tmp}/combined.vcf.gz && /usr/bin/tabix -p vcf ${tmp}/combined.vcf.gz && \
    /opt/conda/bin/python /gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/addReadCountsToVcfCRAM3.py -n ${MinReads} -v ${MinVAF} -r ${refFasta} ${tmp}/combined.vcf.gz ${Bam} ${Name} | \
    bgzip -c > ${Name}.combined_tagged.vcf.gz && /usr/bin/tabix -p vcf ${Name}.combined_tagged.vcf.gz
  }
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "10 G"
    job_group: jobGroup
  }
  output {
    File combined_vcf_file = "${Name}.combined_tagged.vcf.gz"
  }

}

task annotate_variants {
  String Vcf
  String refFasta
  String Vepcache
  String Cytobands
  File CustomAnnotationVcf
  File CustomAnnotationIndex
  String CustomAnnotationParameters
  Float? maxAF
  String Name
  String jobGroup
  String? tmp
  String? docker = "mgibio/basespace_chromoseq:v12"
  
  command {
    set -eo pipefail && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf --vcf --fasta ${refFasta} --hgvs --symbol --term SO --per_gene -o ${Name}.annotated.vcf \
    -i ${Vcf} --custom ${Cytobands},cytobands,bed --custom ${CustomAnnotationVcf},${CustomAnnotationParameters} --offline --cache --max_af --dir ${Vepcache} && \
    /opt/htslib/bin/bgzip -c ${Name}.annotated.vcf > ${Name}.annotated.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/ensembl-vep/filter_vep -i ${Name}.annotated.vcf.gz --format vcf -o ${Name}.annotated_filtered.vcf \
    --filter "(MAX_AF < ${default='0.001' maxAF} or not MAX_AF) or MYELOSEQ_TCGA_AC or MYELOSEQ_MDS_AC" && \
    /opt/htslib/bin/bgzip -c ${Name}.annotated_filtered.vcf > ${Name}.annotated_filtered.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz
  }
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File annotated_vcf = "${Name}.annotated.vcf.gz"
    File annotated_filtered_vcf = "${Name}.annotated_filtered.vcf.gz"
  }
}

task annotate_svs {
  String Vcf
  String CNV
  String refFasta
  String refFastaIndex
  String Vepcache
  String Name
  String gender
  String jobGroup
  String SVAnnot
  String Translocations
  String Cytobands
  Int? minCNAsize
  Float? minCNAabund
  Int? lowCNAsize
  Float? lowCNAabund
  
  String? tmp
  String? docker = "mgibio/basespace_chromoseq:v12"
  
  command {
    set -eo pipefail && \
    perl /gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/ichorToVCF.pl -g ${gender} -minsize ${minCNAsize} \
    -minabund ${minCNAabund} -lowsize ${lowCNAsize} \
    -lowabund ${lowCNAabund} -r ${refFasta} ${CNV} | bgzip -c > cnv.vcf.gz && \
    /opt/htslib/bin/tabix -p vcf cnv.vcf.gz && \
    /opt/conda/envs/python2/bin/bcftools query -l cnv.vcf.gz > name.txt && \
    perl /gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/FilterManta.pl -a ${minCNAabund} -r ${refFasta} -k ${Translocations} ${Vcf} filtered.vcf && \
    /opt/conda/envs/python2/bin/svtools afreq filtered.vcf | \
    /opt/conda/envs/python2/bin/svtools vcftobedpe -i stdin | \
    /opt/conda/envs/python2/bin/svtools varlookup -d 200 -c BLACKLIST -a stdin -b ${SVAnnot} | \
    /opt/conda/envs/python2/bin/svtools bedpetovcf | \
    /usr/local/bin/bedtools sort -header -g ${refFastaIndex} -i stdin | bgzip -c > filtered.tagged.vcf.gz && \
    /opt/conda/envs/python2/bin/bcftools reheader -s name.txt filtered.tagged.vcf.gz > filtered.tagged.reheader.vcf.gz && \
    /opt/htslib/bin/tabix -p vcf filtered.tagged.reheader.vcf.gz && \
    /opt/conda/envs/python2/bin/bcftools concat -a cnv.vcf.gz filtered.tagged.reheader.vcf.gz | \
    /usr/local/bin/bedtools sort -header -g ${refFastaIndex} -i stdin > svs.vcf && \
    /opt/conda/envs/python2/bin/python /usr/local/src/manta/libexec/convertInversion.py /usr/local/bin/samtools ${refFasta} svs.vcf | bgzip -c > svs.vcf.gz && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl --format vcf --vcf --fasta ${refFasta} --per_gene --symbol --term SO -o ${Name}.svs_annotated.vcf -i svs.vcf.gz --custom ${Cytobands},cytobands,bed --offline --cache --dir ${Vepcache} && \
    /opt/htslib/bin/bgzip -c ${Name}.svs_annotated.vcf > ${Name}.svs_annotated.vcf.gz && \
    /opt/htslib/bin/tabix -p vcf ${Name}.svs_annotated.vcf.gz
  }
  
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "24 G"
    job_group: jobGroup
  }
  
  output {
    File vcf = "${Name}.svs_annotated.vcf.gz"
    File vcf_index = "${Name}.svs_annotated.vcf.gz.tbi"
  }
}


task make_report {
  String SVVCF
  String GeneVCF
  String KnownGenes
  String MappingSummary
  String? CoverageSummary
  String SVQC
  String GeneQC
  String Name
  String jobGroup
  String tmp
  String? docker = ""
  Int? MinGeneCov
  Int? MinFracGene20
  Int? MinRegionCov
  Int? MinFracRegion10
  
  command <<<
    cat ${MappingSummary} ${CoverageSummary} | cut -d ',' -f 3,4 | sort -u > qc.txt && \
    /opt/conda/bin/python /gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/make_report3.py ${Name} ${GeneVCF} ${SVVCF} ${KnownGenes} "qc.txt" ${GeneQC} ${SVQC} > "${Name}.chromoseq.txt"
  >>>
  
  runtime {
    docker_image: docker
    job_group: jobGroup
  }
  
  output {
    File report = "${Name}.chromoseq.txt"
  }
  
}

task make_igv {
  String Name
  String docker
  
  command {
    cat <<EOF > ${Name}.igv.xml
    <?xml version="1.0" encoding="UTF-8"?>
    <Session genome="hg38" locus="All" version="3">
    <Resources>
    <Resource name="Structural variants" path="${Name}.svs_annotated.vcf.gz"/>
    <Resource name="Gene variants" path="${Name}.annotated_filtered.vcf.gz"/>
    <Resource name="Log2Ratio CN" path="${Name}.l2r.bw"/>
    <Resource name="Copy Number Est." path="${Name}.cn.bw"/>
    <Resource name="Copy Number Call" path="${Name}.cnv.bed"/>
    <Resource name="Ensemble Genes" path="http://www.broadinstitute.org/igvdata/annotations/hg38/EnsemblGenes.ensGene"/>
    </Resources>
    </Session>
    EOF    
  }
  
  runtime {
    docker_image: docker
  }
  
  output {
    File igv_xml = "${Name}.igv.xml"
  }  
}

task remove_files {
  Array[String] files
  String order_by
  String jobGroup
  String docker
  
  command {
    /bin/rm ${sep=" " files}
  }
  runtime {
    docker_image: docker
    job_group: jobGroup
  }
  output {
    String done = stdout()
  }
}

task gather_files {
  Array[String] OutputFiles
  String OutputDir
  String jobGroup
  String? docker = "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
  
  command {
    /bin/mv -f -t ${OutputDir}/ ${sep=" " OutputFiles}
  }
  runtime {
    docker_image: "ubuntu:xenial"
  }
  output {
    String done = stdout()
  }
}

task return_object {
  Array[Object] obj
  command {
    cat ${write_objects(obj)} > "obj.tsv"
  }
  
  output {
    File results = "obj.tsv"
  }
  
}
