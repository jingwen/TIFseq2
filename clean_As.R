library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
args <- commandArgs(TRUE)
# clean the internal A riched regions in human genome
clean_Arich <- function(file_name,extend_size, maxA){
  reads=read.table(file=file_name,sep="\t",col.names = c("seqname","pos","strand","score"))
  polyA_site <- GRanges(seqnames = reads$seqname,ranges=IRanges(reads$pos,reads$pos),strand=reads$strand, score=reads$score)
  cat(file_name,"cluster before filtering:",length(polyA_site),"\n")
  names(polyA_site)=paste(reads$seqname,reads$pos,reads$strand,sep=":")
  polyA_site=flank(polyA_site,extend_size,start=FALSE)
  genome=BSgenome.Hsapiens.UCSC.hg38
  v=Views(genome,polyA_site)
  polyA_freq=alphabetFrequency(v,baseOnly=TRUE)
  rownames(polyA_freq)=names(polyA_site)
  
  polyA_site=polyA_site[which(polyA_freq[,1] <= maxA),]
  read=str_split_fixed(names(polyA_site),":",3)
  polyA_site=flank(polyA_site,1)
  cat("cluster after filtering:",length(polyA_site),"\n")
  df <- data.frame(seqnames=seqnames(polyA_site),starts=start(polyA_site),strands=strand(polyA_site),scores=polyA_site$score)
  write.table(df,file=file_name,quote=F, sep="\t", row.names=F, col.names=F)
}
dir_3end=args[1]
for (files in list.files(dir_3end,full.names=TRUE)){
  clean_Arich(files,10,6)
}
