load("cluster.rda")
library(data.table)
library(GenomicRanges)

## Preparation of read pair count
addCountColumn <- function(Count.all.samples, count) {
  if(is.null(Count.all.samples))
    return(count)
  merge(Count.all.samples, count, all.x = TRUE, all.y = TRUE)
}

readCount <- function(count_files){
  sampleLabels=sub("_boundary.txt","",basename(count_files))
  Count.all.samples=NULL
  for (i in 1:length(count_files)){
    count=read.table(file=count_files[i], header=TRUE, col.names=c("chr", "end5", "end3", "strand", sampleLabels[i]),
                     colClasses = c("character", "integer", "integer","character", "integer"))
    count=data.table(count)
    setkeyv(count,cols=c("chr","end5","end3","strand"))
    Count.all.samples=addCountColumn(Count.all.samples,count)
  }
  Count.all.samples[is.na(Count.all.samples)]<-0
  return(Count.all.samples)
}

## assign read pairs into cluster pairs
assign2cluster <- function(end5,end3, Count.all.samples){
  TSS=CTSScoordinatesGR(end5)
  TSS=TSS[TSS$cluster!=""]
  TES=CTSScoordinatesGR(end3)
  TES=TES[TES$cluster!=""]
  boundary5=data.frame(chr=TSS@seqnames,start=TSS@ranges,strand=TSS@strand,cluster=TSS@elementMetadata$cluster)
  boundary5=data.table(boundary5)
  colnames(boundary5)=c("chr","end5","strand","cluster5")
  setkeyv(boundary5,cols=c("chr","end5","strand"))
  boundary3=data.frame(chr=TES@seqnames,end=TES@ranges,strand=TES@strand,cluster=TES@elementMetadata$cluster)
  boundary3=data.table(boundary3)
  colnames(boundary3)=c("chr","end3","strand","cluster3")
  setkeyv(boundary3,cols=c("chr","end3","strand"))
  setkeyv(Count.all.samples,cols=c("chr","end5","strand"))
  Count.all.samples=merge(Count.all.samples,boundary5)
  setkeyv(Count.all.samples,cols=c("chr","end3","strand"))
  Count.all.samples=merge(Count.all.samples,boundary3)
  Count.all.samples=as.data.frame(Count.all.samples)
  return(Count.all.samples)
}
count_files=list.files("~/TIFSeq/K562/STAR_i200K_10ntExon/count_2nd", full.names = TRUE)
sample_size=length(count_files)
Counts <<- readCount(count_files)
Count.samples <<- assign2cluster(K562_end5_dist10,K562_end3_dist10,Counts)

# set the boundaries of isoforms
defineBoundary <- function(count,min_mate,max_mate,min_count=0,peak5,peak3){
  col_end=5+sample_size-1
  cluster_boundary=aggregate(x=count[colnames(count)[5:col_end]],by=count[c("chr","cluster5","cluster3","strand")],sum)
  cluster_boundary$cluster5=as.character(cluster_boundary$cluster5)
  cluster_boundary$cluster3=as.character(cluster_boundary$cluster3)
  cluster_boundary$totalCov=rowSums(cluster_boundary[5:col_end])
  cluster_boundary$peak_5end=peak5[cluster_boundary$cluster5]
  cluster_boundary$peak_3end=peak3[cluster_boundary$cluster3]
  cluster_boundary$start=with(cluster_boundary,ifelse(strand=="+",peak_5end,peak_3end))
  cluster_boundary$end=with(cluster_boundary,ifelse(strand=="+",peak_3end,peak_5end))
  cluster_keep=with(cluster_boundary,cluster_boundary[which(end-start <= max_mate & end-start >= min_mate & totalCov >= min_count),])
  cluster_keep=makeGRangesFromDataFrame(cluster_keep, keep.extra.columns=TRUE)
  names(cluster_keep)=as.factor(cluster_keep)
  return(cluster_keep)
}

library(GenomicRanges)
cluster_pairs=defineBoundary(Count.samples,max_mate=2000000,min_mate=300,min_count=4,peak5=peak_5end_dist10,peak3=peak_3end_dist10)
save(cluster_pairs,file="TIFseq_clusters.rda")
