workflow ChromoSeq {

  String Cram
  String CramIndex
  String Name
  String OutputDir

  String Translocations = "/opt/files/ChromoSeq.translocations.fixed.v2.sorted.hg38.bedpe"
  String SVBed = "/opt/files/ChromoSeq.translocations.qc.bed"
  String Cytobands = "/opt/files/ChromoSeq.hg38.bed"
  String CoverageBed = "/opt/files/GeneCoverageRegions.bed"
  String MantaConfig = "/opt/files/configManta.hg38.py.ini"
  String SVDB = "/opt/files/B38.callset.public.bedpe.gz"

  String Blacklist = "/opt/files/hg38.blacklist.merged.bed"
    
  String Reference = "/gscmnt/gc2555/spencer/refdata/hg38/all_sequences.fa"
  String ReferenceIndex = "/gscmnt/gc2555/spencer/refdata/hg38/all_sequences.fa.fai"
  String ReferenceBED = "/gscmnt/gc2555/spencer/refdata/hg38/all_sequences.fa.bed.gz"
  String Dictionary = "/gscmnt/gc2555/spencer/refdata/hg38/all_sequences.dict"
  String VEP = "/gscmnt/gc2709/info/production_reference_GRCh38DH/CLE/IDTExome/VEP_cache/"

  Float minVarFreq=0.05
  
  String JobGroup = "/dspencer/chromoseq"

  call subset_cram {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Bed=ReferenceBED,
    refFasta=Reference,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call cov_qc as gene_qc {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Name=Name,
    Bed=CoverageBed,
    refFasta=Reference,
    jobGroup=JobGroup
  }

  call cov_qc as sv_qc {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Name=Name,
    Bed=SVBed,
    refFasta=Reference,
    jobGroup=JobGroup
  }
  
  call run_manta {
    input: Bam=Cram,
    BamIndex=CramIndex,
    Config=MantaConfig,
    Reference=Reference,
    ReferenceBED=ReferenceBED,
    SVAnnot=SVDB,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call run_ichor {
    input: Bam=Cram,
    BamIndex=CramIndex,
    refFasta=Reference,
    refIndex=ReferenceIndex,
    ReferenceBED=ReferenceBED,
    Bed=Cytobands,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call run_varscan {
    input: Bam=subset_cram.bamfile,
    BamIndex=subset_cram.bamindex,
    CoverageBed=CoverageBed,
    MinFreq=minVarFreq,
    refFasta=Reference,
    Name=Name,
    jobGroup=JobGroup
  }
  
#  call run_platypus {
#    input: Bam=subset_cram.bamfile,
#    BamIndex=subset_cram.bamindex,
#    CoverageBed=CoverageBed,
#    MinFreq=minVarFreq,
#    Name=Name,
#    refFasta=Reference,
#    jobGroup=JobGroup
#  }
  
  call run_pindel_region as run_pindel_flt3itd {
    input: Bam=subset_cram.bamfile,
    BamIndex=subset_cram.bamindex,
    Reg='chr13:28033987-28034316',
    refFasta=Reference,
    Name=Name,
    jobGroup=JobGroup
  }
  
#  call make_bw {
#    input: in=subset_cram.bamfile,
#    index=subset_cram.bamindex,
#    label=Name,
#    Blacklist=Blacklist,
#    jobGroup=JobGroup
#  }
  
  call combine_variants {
    input: VarscanSNV=run_varscan.varscan_snv_file,
    VarscanIndel=run_varscan.varscan_indel_file,
    PindelITD=run_pindel_flt3itd.pindel_vcf_file,
#    Platypus=run_platypus.platypus_vcf_file,
    Bam=subset_cram.bamfile,
    BamIndex=subset_cram.bamindex,
    refFasta=Reference,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call annotate_variants {
    input: Vcf=combine_variants.combined_vcf_file,
    refFasta=Reference,
    Vepcache=VEP,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call annotate_svs {
    input: Vcf=run_manta.filtered_vcf,
    refFasta=Reference,
    Vepcache=VEP,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call make_report {
    input: VCF=annotate_svs.vcf,
    CNV=run_ichor.report,
    VARS=annotate_variants.annotated_filtered_tsv,
    TranslocationsBED=Translocations,
    CytobandsBED=Cytobands,
    Name=Name,
    jobGroup=JobGroup
  }
  
  call make_igv {
    input: Name=Name
  }
  
  call gather_files {
    input: OutputFiles=[annotate_svs.vcf,
    annotate_svs.vcf_index,
    run_ichor.params,
    run_ichor.seg,
    run_ichor.genomewide_pdf,
    run_ichor.allgenomewide_pdf,
    run_ichor.report,run_ichor.rdata,run_ichor.wig,
    run_ichor.cn_bw,run_ichor.l2r_bw,run_ichor.bed,
    run_ichor.correct_pdf,
    gene_qc.qc_out,
    gene_qc.region_dist,
    gene_qc.global_dist,
    sv_qc.qc_out,
    sv_qc.region_dist,
    annotate_variants.annotated_filtered_vcf,
    annotate_variants.annotated_filtered_tsv,
    make_report.report,  #make_bw.bigwig_file,
    make_igv.igv_xml],
    OutputDir=OutputDir,
    jobGroup=JobGroup
  }
  
  call remove_files as cleanup {
    input: files=[subset_cram.bamfile,subset_cram.bamindex],
    order_by=gather_files.done,
    jobGroup=JobGroup    
  }
  
}

task cov_qc {
  String Cram
  String CramIndex
  String Bed
  String Name
  String refFasta
  String jobGroup
  
  command <<<
    /opt/conda/bin/mosdepth -n -f ${refFasta} -t 4 -i 2 -Q 20 -b ${Bed} --thresholds 10,20,30,40 "${Name}" ${Cram} && \
    /usr/local/bin/bedtools intersect -header -b "${Name}.regions.bed.gz" -a "${Name}.thresholds.bed.gz" -wo | \
    awk -v OFS="\t" '{ if (NR==1){ print $0,"%"$5,"%"$6,"%"$7,"%"$8,"MeanCov"; } else { print $1,$2,$3,$4,$5,$6,$7,$8,sprintf("%.2f\t%.2f\t%.2f\t%.2f",$5/$NF*100,$6/$NF*100,$7/$NF*100,$8/$NF*100),$(NF-1); } }' > "${Name}."$(basename ${Bed} .bed)".covqc.txt" && \
    mv "${Name}.mosdepth.region.dist.txt" "${Name}.mosdepth."$(basename ${Bed} .bed)".region.dist.txt"
  >>>
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
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
  String SVAnnot
  String jobGroup
  
  command <<<
    /usr/local/src/manta/bin/configManta.py --config=/opt/files/configManta.hg38.py.ini --tumorBam=${Bam} --referenceFasta=${Reference} \
    --runDir=manta --callRegions=${ReferenceBED} --outputContig && \
    ./manta/runWorkflow.py -m local -q research-hpc -j 4 -g 32 && \
    /opt/conda/envs/python2/bin/python /usr/local/src/manta/libexec/convertInversion.py /usr/local/bin/samtools ${Reference} ./manta/results/variants/tumorSV.vcf.gz | \
    /opt/conda/envs/python2/bin/svtools afreq | /opt/conda/envs/python2/bin/svtools vcftobedpe -i stdin | \
    /opt/conda/envs/python2/bin/svtools varlookup -d 200 -c POPFREQ -a stdin -b ${SVAnnot} | \
    /opt/conda/envs/python2/bin/svtools bedpetovcf | /opt/conda/envs/python2/bin/svtools vcfsort > ${Name}.tumorSV.vcf && \
    perl /usr/local/bin/BlatContigs.pl -r ${Reference} ${Name}.tumorSV.vcf ${Name}.tumorSV.filtered.vcf && \
    bgzip ${Name}.tumorSV.filtered.vcf && tabix -p vcf ${Name}.tumorSV.filtered.vcf.gz
  >>>
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "4"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File filtered_vcf = "${Name}.tumorSV.filtered.vcf.gz"
    File filtered_index = "${Name}.tumorSV.filtered.vcf.gz.tbi"
  }
}


task run_ichor {
  String Bam
  String BamIndex
  String ReferenceBED
  String Bed
  String refFasta
  String refIndex
  String Name
  String jobGroup

  command <<<
    /usr/local/bin/bedtools makewindows -b ${ReferenceBED} -w 500000 > /tmp/windows.bed && \
    /usr/local/bin/samtools view -f 0x2 -F 0x400 -q 20 -T ${refFasta} -b ${Bam} | \
    /usr/local/bin/intersectBed -sorted -nobuf -c -bed -b stdin -a /tmp/windows.bed | \
    awk -v window=500000 'BEGIN { chr=""; } { if ($1!=chr){ printf("fixedStep chrom=%s start=1 step=%d span=%d\n",$1,window,window); chr=$1; } print $4; }' > "${Name}.tumor.wig" && \
    /usr/local/bin/Rscript /usr/local/bin/runIchorCNA.R \
    --id ${Name} \
    --WIG "${Name}.tumor.wig" --ploidy "c(2)" --normal "c(0.1,0.5,.85)" --maxCN 3 \
    --gcWig /usr/local/lib/R/site-library/ichorCNA/extdata/gc_hg38_500kb.wig \
    --mapWig /usr/local/lib/R/site-library/ichorCNA/extdata/map_hg38_500kb.wig \
    --centromere /usr/local/lib/R/site-library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
    --normalPanel /opt/files/nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.rds \
    --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --fracReadsInChrYForMale 0.0005 \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --txnE 0.999999 --txnStrength 1000000 --genomeStyle UCSC --outDir ./ && \
    awk -v OFS="\t" '$7!=2 && NR>1 { print $2,$3,$4,$5,$6,$7,$8,$9; }' "${Name}.seg.txt" > results.bed && \
    if [[ -s results.bed ]]; then /usr/local/bin/intersectBed -a results.bed -b ${Bed} -wa -wb | \
    /usr/local/bin/bedtools groupby -g 1,2,3,4,5,6,7 -c 12,12,13,13,14 -o distinct,count_distinct,distinct,count_distinct,distinct >> "${Name}.cnv_report.txt"; else touch "${Name}.cnv_report.txt"; fi && \
    awk -v OFS="\t" '{ color="255,0,0"; if ($5>0){ color="0,0,255"; } n=split($8,a,","); print $1,$2,$3,$7"("a[0]"-"a[n]")",".",".",$2,$3,color; }' "${Name}.cnv_report.txt" > "${Name}.cnv.bed" && \
    awk -v OFS="\t" 'NR > 1 && !/NA/ { print $1,$2,$3,$4; }' "${Name}.cna.seg" | sort -k 1,1 -k 2,2n > /tmp/cn.bedgraph && \
    /usr/bin/bedGraphToBigWig /tmp/cn.bedgraph ${refIndex} "${Name}.cn.bw" && \
    awk -v OFS="\t" 'NR > 1 && !/NA/ { print $1,$2,$3,$6; }' "${Name}.cna.seg" | sort -k 1,1 -k 2,2n > /tmp/l2r.bedgraph && \
    /usr/bin/bedGraphToBigWig /tmp/l2r.bedgraph ${refIndex} "${Name}.l2r.bw" && \
    mv ${Name}/*.pdf .
  >>>
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "16 G"
    job_group: jobGroup
  }
  
  output {
    File params = "${Name}.params.txt"
    File seg = "${Name}.seg.txt"
    File cn_bw = "${Name}.cn.bw"
    File l2r_bw = "${Name}.l2r.bw"
    File bed = "${Name}.cnv.bed"       
    File report = "${Name}.cnv_report.txt"
    File genomewide_pdf = "${Name}_genomeWide.pdf"
    File allgenomewide_pdf = "${Name}_genomeWide_all_sols.pdf"
    File correct_pdf = "${Name}_correct.pdf"
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
  String CoverageBed
  String refFasta
  String Name
  String jobGroup
    
  command <<<
    /usr/local/bin/samtools mpileup -f ${refFasta}".gz" -l ${CoverageBed} ${Bam} > /tmp/mpileup.out && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2snp /tmp/mpileup.out --min-coverage ${default=8 MinCov} --min-reads2 ${default=5 MinReads} \
    --min-var-freq ${default="0.02" MinFreq} --output-vcf > ${Name}.snv.vcf && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2indel /tmp/mpileup.out --min-coverage ${default=8 MinCov} --min-reads2 ${default=5 MinReads} \
    --min-var-freq ${default="0.02" MinFreq} --output-vcf > ${Name}.indel.vcf
  >>>
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "2"
    memory: "16 G"
    job_group: jobGroup
  }
  output {
    File varscan_snv_file = "${Name}.snv.vcf"
    File varscan_indel_file = "${Name}.indel.vcf"
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
  
  command <<<
    (set -eo pipefail && /usr/local/bin/samtools view -T ${refFasta}".gz" ${Bam} ${Reg} | /opt/pindel-0.2.5b8/sam2pindel - /tmp/in.pindel ${default=250 Isize} tumor 0 Illumina-PairEnd) && \
    /usr/local/bin/pindel -f ${refFasta} -p /tmp/in.pindel -c ${Reg} -o /tmp/out.pindel && \
    /usr/local/bin/pindel2vcf -P /tmp/out.pindel -G -r ${refFasta} -e ${default=3 MinReads} -R hg38 -d hg38 -v pindel.vcf && \
    /bin/sed 's/END=[0-9]*\;//' pindel.vcf > ${Name}.pindel.vcf
  >>>
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "16 G"
    job_group: jobGroup
  }
  output {
    File pindel_vcf_file = "${Name}.pindel.vcf"
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
  
  command <<<
    /usr/bin/awk '{ print $1":"$2+1"-"$3; }' ${CoverageBed} > "regions.txt" && \
    /opt/conda/bin/octopus -R ${refFasta} -I ${Bam} -t regions.txt -C cancer > "${Name}.vcf" && \
    /bin/sed 's/VCFv4.3/VCFv4.1/' "${Name}.vcf" > "${Name}.platypus.vcf"     
  >>>
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File platypus_vcf_file = "${Name}.platypus.vcf"
  }
}

task subset_cram {
  String Cram
  String CramIndex
  String refFasta
  String Bed
  String Name
  String jobGroup
  
  command {
    /usr/local/bin/samtools view -T ${refFasta} -L ${Bed} -b -o "${Name}.subset.bam" ${Cram} && \
    /usr/local/bin/samtools index "${Name}.subset.bam"
  }
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "16 G"
    job_group: jobGroup
  }
  
  output {
    File bamfile = "${Name}.subset.bam"
    File bamindex = "${Name}.subset.bam.bai"
  }
  
}

task make_bw {
  String in
  String index
  String label
  String Blacklist
  String jobGroup
  
  Int? genome_size
  
  command {
    export PYTHONPATH=/opt/conda/lib/python3.6/site-packages/ && \
    /opt/conda/bin/bamCoverage --bam ${in} -o "${label}.bw" --effectiveGenomeSize ${default=2451960000 genome_size} --normalizeUsing RPGC \
    --ignoreDuplicates -bl ${Blacklist} --binSize 50 --minMappingQuality 1 --extendReads -p 4 -ignore X Y MT
  }
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "4"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File bigwig_file = "${label}.bw"
  }
}           
  

task combine_variants {
  String VarscanSNV
  String VarscanIndel
  String PindelITD
#  String Platypus
  String Bam
  String BamIndex
  String refFasta
  String Name
  String jobGroup
  
  command {
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants -R ${refFasta} --variant:varscanIndel ${VarscanIndel} \
    --variant:varscanSNV ${VarscanSNV} --variant:PindelITD ${PindelITD} -o /tmp/out.vcf --genotypemergeoption UNIQUIFY && \
    /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ${refFasta} --variant /tmp/out.vcf -o /tmp/combined.vcf && \
    /opt/conda/bin/python /usr/local/bin/addReadCountsToVcfCRAM.py -r ${refFasta} /tmp/combined.vcf ${Bam} ${Name} > ${Name}.combined_tagged.vcf
  }
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "10 G"
    job_group: jobGroup
  }
  output {
    File combined_vcf_file = "${Name}.combined_tagged.vcf"
  }
  
}

task annotate_variants {
  String Vcf
  String refFasta
  String Vepcache
  Float? maxAF
  String Name
  String jobGroup
  
  command {
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf --vcf --plugin Downstream --plugin Wildtype --fasta ${refFasta} --hgvs --symbol --term SO --flag_pick -o ${Name}.annotated.vcf \
    -i ${Vcf} --offline --cache --af_gnomad --dir ${Vepcache} && \
    /opt/htslib/bin/bgzip -c ${Name}.annotated.vcf > ${Name}.annotated.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/ensembl-vep/filter_vep -i ${Name}.annotated.vcf.gz --format vcf -o ${Name}.annotated_filtered.vcf \
    --filter "(gnomAD_AF < 0.001 < 0.001 and gnomAD_AFR_AF < 0.001 and gnomAD_SAS_AF < 0.001 and gnomAD_EAS_AF < 0.001 and gnomAD_NFE_AF < 0.001 and gnomAD_AMR_AF < 0.001 and gnomAD_OTH_AF < 0.001 and gnomAD_FIN_AF < 0.001) or not gnomAD_AF" && \
    /opt/htslib/bin/bgzip -c ${Name}.annotated_filtered.vcf > ${Name}.annotated_filtered.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz && \
    /usr/bin/java -Xmx4g -jar /opt/GenomeAnalysisTK.jar -T VariantsToTable \
    -R ${refFasta} --variant ${Name}.annotated_filtered.vcf.gz -o ${Name}.variants.tsv \
    -F CHROM -F POS -F ID -F REF -F ALT -F set \
    -GF GT -GF NR -GF NV -GF VAF && \
    /opt/conda/bin/python /usr/local/bin/add_annotations_to_table_helper.py ${Name}.variants.tsv ${Name}.annotated_filtered.vcf.gz Consequence,SYMBOL,Feature_type,Feature,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,HGNC_ID,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,CLIN_SIG,SOMATIC,PHENO ./ && \
    mv variants.annotated.tsv ${Name}.variants_annotated.tsv
    
  }
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "32 G"
    job_group: jobGroup
  }
  output {
    File annotated_vcf = "${Name}.annotated.vcf.gz"
    File annotated_filtered_vcf = "${Name}.annotated_filtered.vcf.gz"
    File annotated_filtered_tsv = "${Name}.variants_annotated.tsv"
  }
}

task annotate_svs {
  String Vcf
  String refFasta
  String Vepcache
  String Name
  String jobGroup
  
  command {
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf --vcf --plugin Downstream --plugin Wildtype --fasta ${refFasta} --symbol --term SO --flag_pick -o ${Name}.svs_annotated.vcf \
    -i ${Vcf} --offline --cache --dir ${Vepcache} && \
    /opt/htslib/bin/bgzip -c ${Name}.svs_annotated.vcf > ${Name}.svs_annotated.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.svs_annotated.vcf.gz
  }
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    cpu: "1"
    memory: "10 G"
    job_group: jobGroup
  }
  
  output {
    File vcf = "${Name}.svs_annotated.vcf.gz"
    File vcf_index = "${Name}.svs_annotated.vcf.gz.tbi"
  }
}


task make_report {
  String VCF
  String CNV
  String VARS
  String TranslocationsBED
  String CytobandsBED
  String Name
  String jobGroup
  
  command {
    perl /usr/local/bin/ChromoSeqReporter.hg38.pl ${Name} ${VARS} ${CNV} ${VCF} > "${Name}.chromoseq.txt"
  }
  
  runtime {
    docker_image: "dhspence/docker-chromoseq"
    job_group: jobGroup
  }
  
  output {
    File report = "${Name}.chromoseq.txt"
  }
  
}

task make_igv {
  String Name
  
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
    docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
  }
  
  output {
    File igv_xml = "${Name}.igv.xml"
  }  
}

task remove_files {
  Array[String] files
  String order_by
  String jobGroup
  
  command {
    /bin/rm ${sep=" " files}
  }
  runtime {
    docker_image: "ubuntu:xenial"
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
