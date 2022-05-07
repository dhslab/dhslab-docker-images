workflow ChromoSeq {

  # required inputs
  
  String Cram 
  String CramIndex  
  String Name  
  String OutputDir 

  # the rest are defaults
  
  String Translocations = "/opt/files/chromoseq_translocations.bedpe"
  String GenesBed = "/opt/files/chromoseq_genes.bed"
  
  String Cytobands = "/opt/files/hg38.cytoBandIdeo.bed.gz"
  String SVDB = "/opt/files/chromoseq_sv_filter.bedpe.gz"
  
  String MantaConfig = "/opt/files/configManta.hg38.py.ini"
  
  String Reference = "/opt/files/all_sequences.fa"
  String ReferenceIndex = "/opt/files/all_sequences.fa.fai"
  String ReferenceBED = "/opt/files/all_sequences.fa.bed.gz"
  String VEP = "/opt/files"

  String gcWig = "/usr/local/lib/R/site-library/ichorCNA/extdata/gc_hg38_500kb.wig"
  String mapWig = "/usr/local/lib/R/site-library/ichorCNA/extdata/map_hg38_500kb.wig"
  String ponRds = "/opt/files/nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds"
  String centromeres = "/usr/local/lib/R/site-library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt"  
  
  String tmp = "/tmp"

  String GenomeBuild = "hg38"
  
  Float minVarFreq=0.02
  
  call prepare_bed {
    input: Bedpe=Translocations,
    Bed=GenesBed,
    Reference=ReferenceBED,
    tmp=tmp
  }
    
  call cov_qc as gene_qc {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Name=Name,
    Bed=GenesBed,
    refFasta=Reference,
    tmp=tmp
  }

  call cov_qc as sv_qc {
    input: Cram=Cram,
    CramIndex=CramIndex,
    Name=Name,
    Bed=prepare_bed.svbed,
    refFasta=Reference,
    tmp=tmp
  }
  
  call run_manta {
    input: Bam=Cram,
    BamIndex=CramIndex,
    Config=MantaConfig,
    Reference=Reference,
    ReferenceBED=ReferenceBED,
    Name=Name,
    tmp=tmp
  }
    
  scatter (chr in prepare_bed.chroms){
    call count_reads {
      input: Bam=Cram,
      BamIndex=CramIndex,
      ReferenceBED=ReferenceBED,
      refFasta=Reference,
      refIndex=ReferenceIndex,
      Chrom=chr,
      tmp=tmp
    }
  }
  
  call run_ichor {
    input: Bam=Cram,
    BamIndex=CramIndex,
    refFasta=Reference,
    ReferenceBED=ReferenceBED,
    CountFiles=count_reads.counts_bed,
    gcWig=gcWig,
    mapWig=mapWig,
    ponRds=ponRds,
    centromeres=centromeres,
    Name=Name,
    tmp=tmp
  }
  
  call run_varscan {
    input: Bam=Cram,
    BamIndex=CramIndex,
    CoverageBed=GenesBed,
    MinFreq=minVarFreq,
    refFasta=Reference,
    Name=Name,
    tmp=tmp
  }
  
  call run_pindel_region as run_pindel_flt3itd {
    input: Bam=Cram,
    BamIndex=CramIndex,
    Reg='chr13:28033987-28034316',
    refFasta=Reference,
    genome=GenomeBuild,
    Name=Name,
    tmp=tmp
  }
    
  call combine_variants {
    input: VCFs=[run_varscan.varscan_snv_file,
    run_varscan.varscan_indel_file,
    run_pindel_flt3itd.pindel_vcf_file],
    Bam=Cram,
    BamIndex=CramIndex,
    refFasta=Reference,
    Name=Name,
    tmp=tmp
  }
  
  call annotate_variants {
    input: Vcf=combine_variants.combined_vcf_file,
    refFasta=Reference,
    Vepcache=VEP,
    Cytobands=Cytobands,
    Name=Name,
    tmp=tmp
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
    Name=Name,
    tmp=tmp
  }
  
  call make_report {
    input: SVVCF=annotate_svs.vcf,
    GeneVCF=annotate_variants.annotated_filtered_vcf,
    KnownGenes=prepare_bed.genes,
    GeneQC=gene_qc.qc_out,
    SVQC=sv_qc.qc_out,
    CovQC=gene_qc.global_dist,
    Name=Name,
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
    OutputDir=OutputDir
  }
  
}

task prepare_bed {
  String Bedpe
  String Bed
  String Reference
  String? tmp
  
  command <<<
    awk -v OFS="\t" '{ split($7,a,"_"); print $1,$2,$3,a[1],".",$9; print $4,$5,$6,a[2],".",$10; }' ${Bedpe} | sort -u -k 1,1V -k 2,2n > sv.bed
    cat sv.bed ${Bed} | cut -f 4 > genes.txt
    gunzip -c ${Reference} | cut -f 1 > chroms.txt
  >>>

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
  String tmp
  
  command <<<
    set -eo pipefail && \
    /opt/conda/bin/mosdepth -n -f ${refFasta} -t 4 -i 2 -x -Q 20 -b ${Bed} --thresholds 10,20,30,40 "${Name}" ${Cram} && \
    /usr/local/bin/bedtools intersect -header -b "${Name}.regions.bed.gz" -a "${Name}.thresholds.bed.gz" -wo | \
    awk -v OFS="\t" '{ if (NR==1){ print $0,"%"$5,"%"$6,"%"$7,"%"$8,"MeanCov"; } else { print $1,$2,$3,$4,$5,$6,$7,$8,sprintf("%.2f\t%.2f\t%.2f\t%.2f",$5/$NF*100,$6/$NF*100,$7/$NF*100,$8/$NF*100),$(NF-1); } }' > "${Name}."$(basename ${Bed} .bed)".covqc.txt" && \
    mv "${Name}.mosdepth.region.dist.txt" "${Name}.mosdepth."$(basename ${Bed} .bed)".region.dist.txt"
  >>>
  
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
  String tmp
  
  command <<<
    set -eo pipefail && \
    /usr/local/src/manta/bin/configManta.py --config=${Config} --tumorBam=${Bam} --referenceFasta=${Reference} \
    --runDir=manta --callRegions=${ReferenceBED} --outputContig && \
    ./manta/runWorkflow.py -m local -q research-hpc -j 4 -g 32 && \
    zcat ./manta/results/variants/tumorSV.vcf.gz | /bin/sed 's/DUP:TANDEM/DUP/g' > fixed.vcf && \
    /usr/local/bin/duphold_static -v fixed.vcf -b ${Bam} -f ${Reference} -t 4 -o ${Name}.tumorSV.vcf && \
    bgzip ${Name}.tumorSV.vcf && /usr/bin/tabix ${Name}.tumorSV.vcf.gz
  >>>
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
  String refFasta
  String refIndex
  String tmp
  
  command {
    set -eo pipefail && \
    (/usr/local/bin/bedtools makewindows -b ${ReferenceBED} -w 500000 | \
    awk -v OFS="\t" -v C="${Chrom}" '$1==C && NF==3' > ${tmp}/${Chrom}.windows.bed) && \
    /usr/local/bin/samtools view -b -f 0x2 -F 0x400 -q 20 -T ${refFasta} ${Bam} ${Chrom} | \
    /usr/local/bin/intersectBed -sorted -nobuf -c -bed -b stdin -a ${tmp}/${Chrom}.windows.bed > ${Chrom}.counts.bed
  }
  output {
    File counts_bed = "${Chrom}.counts.bed"
  }
}


task run_ichor {
  String Bam
  String BamIndex
  String ReferenceBED
  Array[String] CountFiles
  String refFasta
  String Name
  String genomeStyle = "UCSC"
  String gcWig
  String mapWig
  String ponRds
  String centromeres
  String? tmp
  
  command <<<
    set -eo pipefail && \
    cat ${sep=" " CountFiles} | sort -k 1,1V -k 2,2n | \
    awk -v window=500000 'BEGIN { chr=""; } { if ($1!=chr){ printf("fixedStep chrom=%s start=1 step=%d span=%d\n",$1,window,window); chr=$1; } print $4; }' > "${Name}.tumor.wig" && \
    /usr/local/bin/Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R --id ${Name} \
    --WIG "${Name}.tumor.wig" --ploidy "c(2)" --normal "c(0.1,0.5,.85)" --maxCN 3 \
    --gcWig ${gcWig} \
    --mapWig ${mapWig} \
    --centromere ${centromeres} \
    --normalPanel ${ponRds} \
    --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --fracReadsInChrYForMale 0.0005 \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --txnE 0.999999 --txnStrength 1000000 --genomeStyle ${genomeStyle} --outDir ./ --libdir /usr/local/bin/ichorCNA/ && \
    mv ${Name}/*.pdf .
  >>>
  
  output {
    File params = "${Name}.params.txt"
    File seg = "${Name}.seg.txt"
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
  String? tmp
  
  command <<<
    /usr/local/bin/samtools mpileup -f ${refFasta} -l ${CoverageBed} ${Bam} > ${tmp}/mpileup.out && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2snp ${tmp}/mpileup.out --min-coverage ${default=6 MinCov} --min-reads2 ${default=3 MinReads} \
    --min-var-freq ${default="0.02" MinFreq} --output-vcf | bgzip -c > ${Name}.snv.vcf.gz && tabix ${Name}.snv.vcf.gz && \
    java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2indel ${tmp}/mpileup.out --min-coverage ${default=6 MinCov} --min-reads2 ${default=3 MinReads} \
    --min-var-freq ${default="0.02" MinFreq} --output-vcf | bgzip -c > ${Name}.indel.vcf.gz && tabix ${Name}.indel.vcf.gz
  >>>
  
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
  String? tmp
  String genome
  
  command <<<
    (set -eo pipefail && /usr/local/bin/samtools view -T ${refFasta} ${Bam} ${Reg} | /opt/pindel-0.2.5b8/sam2pindel - ${tmp}/in.pindel ${default=250 Isize} tumor 0 Illumina-PairEnd) && \
    /usr/local/bin/pindel -f ${refFasta} -p ${tmp}/in.pindel -c ${Reg} -o ${tmp}/out.pindel && \
    /usr/local/bin/pindel2vcf -P ${tmp}/out.pindel -G -r ${refFasta} -e ${default=3 MinReads} -R ${default="hg38" genome} -d ${default="hg38" genome} -v ${tmp}/pindel.vcf && \
    /bin/sed 's/END=[0-9]*\;//' ${tmp}/pindel.vcf | bgzip -c > ${Name}.pindel.vcf.gz && tabix ${Name}.pindel.vcf.gz
  >>>
  
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
  String? tmp
  
  command <<<
    /usr/bin/awk '{ print $1":"$2+1"-"$3; }' ${CoverageBed} > "regions.txt" && \
    /opt/conda/bin/octopus -R ${refFasta} -I ${Bam} -t regions.txt -C cancer > "${Name}.vcf" && \
    /bin/sed 's/VCFv4.3/VCFv4.1/' "${Name}.vcf" > "${Name}.platypus.vcf"     
  >>>
  
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
  String? tmp
  
  command {
    /opt/conda/envs/python2/bin/bcftools merge --force-samples -O z ${sep=" " VCFs} | \
    /opt/conda/envs/python2/bin/bcftools norm -f ${refFasta} > ${tmp}/combined.vcf && \
    /opt/conda/bin/python /usr/local/bin/addReadCountsToVcfCRAM.py -r ${refFasta} ${tmp}/combined.vcf ${Bam} ${Name} | \
    bgzip -c > ${Name}.combined_tagged.vcf.gz && /usr/bin/tabix -p vcf ${Name}.combined_tagged.vcf.gz
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
  Float? maxAF
  String Name
  String? tmp
  
  command {
    set -eo pipefail && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf --vcf --fasta ${refFasta} --hgvs --symbol --term SO --per_gene -o ${Name}.annotated.vcf \
    -i ${Vcf} --custom ${Cytobands},cytobands,bed --offline --cache --max_af --dir ${Vepcache} && \
    /opt/htslib/bin/bgzip -c ${Name}.annotated.vcf > ${Name}.annotated.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/ensembl-vep/filter_vep -i ${Name}.annotated.vcf.gz --format vcf -o ${Name}.annotated_filtered.vcf \
    --filter "(MAX_AF < ${default='0.001' maxAF} or not MAX_AF)" && \
    /opt/htslib/bin/bgzip -c ${Name}.annotated_filtered.vcf > ${Name}.annotated_filtered.vcf.gz && \
    /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz
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
  String SVAnnot
  String Translocations
  String Cytobands
  String? tmp
  
  command {
    set -eo pipefail && \
    perl /usr/local/bin/ichorToVCF.pl -r ${refFasta} ${CNV} | bgzip -c > cnv.vcf.gz && \
    /opt/htslib/bin/tabix -p vcf cnv.vcf.gz && \
    /opt/conda/envs/python2/bin/bcftools query -l cnv.vcf.gz > name.txt && \
    perl /usr/local/bin/FilterManta.pl -r ${refFasta} -k ${Translocations} ${Vcf} filtered.vcf && \
    /opt/conda/envs/python2/bin/svtools afreq filtered.vcf | \
    /opt/conda/envs/python2/bin/svtools vcftobedpe -i stdin | \
    /opt/conda/envs/python2/bin/svtools varlookup -d 200 -c BLACKLIST -a stdin -b ${SVAnnot} | \
    /opt/conda/envs/python2/bin/svtools bedpetovcf | \
    /usr/local/bin/bedtools sort -header -g ${refFastaIndex} -i stdin | bgzip -c > filtered.tagged.vcf.gz && \
    /opt/conda/envs/python2/bin/bcftools reheader -s name.txt filtered.tagged.vcf.gz > filtered.tagged.reheader.vcf.gz && \
    /opt/htslib/bin/tabix -p vcf filtered.tagged.reheader.vcf.gz && \
    /opt/conda/envs/python2/bin/bcftools concat -a cnv.vcf.gz filtered.tagged.reheader.vcf.gz | \
    /usr/local/bin/bedtools sort -header -g ${refFastaIndex} -i stdin | bgzip -c > svs.vcf.gz && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl --format vcf --vcf --fasta ${refFasta} --per_gene --symbol --term SO -o ${Name}.svs_annotated.vcf -i svs.vcf.gz --custom ${Cytobands},cytobands,bed --offline --cache --dir ${Vepcache} && \
    /opt/htslib/bin/bgzip -c ${Name}.svs_annotated.vcf > ${Name}.svs_annotated.vcf.gz && \
    /opt/htslib/bin/tabix -p vcf ${Name}.svs_annotated.vcf.gz
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
  String SVQC
  String GeneQC
  String CovQC
  String Name
  String tmp
  Int? MinGeneCov
  Int? MinFracGene20
  Int? MinRegionCov
  Int? MinFracRegion10
  
  command <<<
    /opt/conda/bin/python /usr/local/bin/make_report.py ${Name} ${GeneVCF} ${SVVCF} ${KnownGenes} > "${Name}.chromoseq.txt" && \
    awk -v F=${default=95 MinFracGene20} -v C=${default=20 MinGeneCov} 'BEGIN { printf("\n*** Gene Coverage Metrics: Exons with <%d%% at 20x or <%dx mean coverage ***\n",F,C) } NR==1 || $10<F || $13<C { print $0; } END { printf("\n"); }' ${GeneQC} >> "${Name}.chromoseq.txt" && \
    awk -v F=${default=95 MinFracRegion10} -v C=${default=10 MinRegionCov} 'BEGIN { printf("\n*** SV Region Coverage Metrics: Genes with <%d%% at 10x or <%dx mean coverage ***\n",F,C) } NR==1 || $9<F || $13<C { print $0; } END { printf("\n"); }' ${SVQC} >> "${Name}.chromoseq.txt"
  >>>
  
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
  
  output {
    File igv_xml = "${Name}.igv.xml"
  }  
}

task remove_files {
  Array[String] files
  String order_by
  
  command {
    /bin/rm ${sep=" " files}
  }
  output {
    String done = stdout()
  }
}

task gather_files {
  Array[String] OutputFiles
  String OutputDir
  
  command {
    /bin/mv -f -t ${OutputDir}/ ${sep=" " OutputFiles}
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
