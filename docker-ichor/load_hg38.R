genomeBuild = "hg38"
genomeStyle = "UCSC"
library(GenomeInfoDb)
bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
  if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
    seqinfo <- Seqinfo(genome=genomeBuild)
  } else {
    seqinfo <- seqinfo(get(bsg))
  }
seqlevelsStyle(seqinfo) <- genomeStyle
seqinfo <- keepSeqlevels(seqinfo, value = paste0("chr",c(1:22,"X")))
saveRDS(seqinfo, file = "/opt/ichorcna/seqinfo_hg38_ucsc.rds")
