#!/bin/bash

echo "start loading modules"
module load bioinfo-tools
module load star

while getopts "hR:I:O:p:j:m:" opt; do
        case $opt in
	h)
		echo "Usage: STAR_align.sh -R <reference genome dir> -I <input fastq dir> -O <output fastq dir> -p <thread number>"
                exit 0
                ;;
        R)genomedir=$OPTARG;;
	I)indir=$OPTARG;;
        O)outdir=$OPTARG;;
        p)thread=$OPTARG;;
	j)intron=$OPTARG;;
	m)mate=$OPTARG;;
	esac
done
echo "Read parameters... DONE"
if [ -z "$thread" ]; then thread=1; fi
if [ -z "$intron" ]; then intron=0; fi
if [ -z "$mate" ]; then mate=0; fi

echo "Referenece: $genomedir"
echo "Input directory: $indir"
echo "Output directory: $outdir"
echo "#thread for alignment: $thread"
echo "maximum intron length: $intron; maximum mate pair distance: $mate"

if [ ! -d "$indir" ];then
        echo "Input directory does not exist!"
        exit 0
fi

if [ ! -d "$outdir" ];then
        echo "Output directory does not exist! Generate directory."
        mkdir $outdir
fi

indir=${indir%/}
outdir=${outdir%/}
for cut5 in $indir/*5cut*.fastq.gz; do
	cut3=${cut5/_5cut/_3cut}
	file=${cut5##*/}
	IFS=_
	name=($file)
	unset IFS
	echo "star --runThreadN $thread --runMode alignReads --genomeDir $genomedir --readFilesIn $cut5 $cut3 --alignIntronMax $intron --alignMatesGapMax $mate --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/${name[0]}."
	star --runThreadN $thread --runMode alignReads --genomeDir $genomedir --readFilesIn $cut5 $cut3 --alignIntronMax $intron --alignMatesGapMax $mate --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/${name[0]}.
done
