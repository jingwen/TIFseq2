library("MultiAssayExperiment")
library("SummarizedExperiment")
library(CAGEr)
library(stringr)
args <- commandArgs(TRUE)
cluster<-function(input_dir,a,r1,r2,n,dis,suffix){
  input_files=list.files(input_dir,full.names=TRUE)
  cluster_end<-CAGEexp(genomeName="BSgenome.Hsapiens.UCSC.hg38",inputFilesType="ctss",
                        inputFiles=input_files, 
                        sampleLabels=sub("-","",sub(suffix,"", basename(input_files))))
  getCTSS(cluster_end)
  normalizeTagCount(cluster_end, method = "powerLaw", 
                    fitInRange = c(r1, r2),alpha = a, T = n)
  clusterCTSS(object=cluster_end,threshold = 1,thresholdIsTpm = TRUE,
              nrPassThreshold = 2, method = "distclu", maxDist = dis,
              removeSingletons = TRUE, keepSingletonsAbove = 1)
  cumulativeCTSSdistribution(cluster_end,clusters = "tagClusters")
  quantilePositions(cluster_end,clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  aggregateTagClusters(cluster_end,tpmThreshold = 1, 
                       excludeSignalBelowThreshold =FALSE, maxDist = 10)
  cumulativeCTSSdistribution(cluster_end,clusters = "consensusClusters")
  quantilePositions(cluster_end,clusters = "consensusClusters", qLow = 0.1, qUp = 0.9)
  return(cluster_end)
}

K562_end5_dist10=cluster(args[1],a=1.49,r1=1,r2=10000,n=4.2*10^6)
K562_end3_dist10=cluster(args[2],a=1.17,r1=1,r2=10000,n=3.8*10^6)
K562_3Tfill_dist10=cluster(args[3],a=1.05,r1=1,r2=100000,n=2.5*10^6,dis=10,suffix="_3end.ctss")

## determined the peak of each consensus cluster
dominant <- function(cluster_end){
  tagNorm=CTSSnormalizedTpm(cluster_end)
  tag=CTSScoordinatesGR(cluster_end)
  tag$score=rowSums(tagNorm[,4:length(tagNorm)])
  tag <- tag[tag$cluster != ""]
  tag_frame=as.data.frame(tag)
  tag_list=split(tag_frame,tag_frame$cluster)
  peak=lapply(tag_list,function(x) x[which.max(x$score),]$pos)
  return(peak)
}

peak_5end_dist10=unlist(dominant(K562_end5_dist10))
peak_3end_dist10=unlist(dominant(K562_end3_dist10))

save(K562_3Tfill_dist10,K562_end5_dist10,K562_end3_dist10,peak_5end_dist10,peak_3end_dist10,file="cluster.rda")
