#!/bin/bash

module load bioinfo-tools
module load star/2.5.3a

while getopts "hR:A:I:O:p:j:m:" opt; do
        case $opt in
	h)
		echo "Usage: STAR_align.sh -R <reference genome dir> -A <annotation gtf> -I <input fastq dir> -O <output fastq dir> -p <thread number> -j <max intron> -m <max mate distance>"
                exit 0
                ;;
        R)genomedir=$OPTARG;;
	A)annotation=$OPTARG;;
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
echo "Annotation file: $annotation"
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

pyenv local anaconda2-2018.12

for cut5 in $indir/*5cut*.fastq.gz; do
	cut3=${cut5/_5cut/_3cut}
	echo "5'end file: $cut5"
	echo "3'end file: $cut3"
	file=${cut5##*/}
	IFS=_
	name=($file)
	echo "${name[0]}"
	unset IFS
#	star --runThreadN $thread --runMode alignReads --genomeDir $genomedir --sjdbGTFfile $annotation --readFilesIn $cut5 $cut3 --alignIntronMax $intron --alignMatesGapMax $mate --alignEndsType Extend5pOfReads12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/${name[0]}.
	star --runThreadN $thread --runMode alignReads --genomeDir $genomedir --sjdbGTFfile $annotation --readFilesIn $cut5 --alignIntronMax $intron --alignEndsType Extend5pOfRead1 --alignSJoverhangMin 10 --outFilterMismatchNmax 5 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMattributes NH HI NM AS --outSAMattrRGline ID:${name[0]}_5end --outFileNamePrefix $outdir/${name[0]}_5end.
	samtools sort -n -o $outdir/${name[0]}_5end.name.bam $outdir/${name[0]}_5end.Aligned.out.bam
	star --runThreadN $thread --runMode alignReads --genomeDir $genomedir --sjdbGTFfile $annotation --readFilesIn $cut3 --alignIntronMax $intron --alignEndsType Extend5pOfRead1 --alignSJoverhangMin 10 --outFilterMismatchNmax 5 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMattributes NH HI NM AS --outSAMattrRGline ID:${name[0]}_3end --outFileNamePrefix $outdir/${name[0]}_3end.
	samtools sort -@ $thread -n -o $outdir/${name[0]}_3end.name.bam $outdir/${name[0]}_3end.Aligned.out.bam
	samtools merge -@ $thread -f -n $outdir/${name[0]}.name.bam $outdir/${name[0]}_5end.name.bam $outdir/${name[0]}_3end.name.bam
	rm $outdir/${name[0]}_*end.name.bam
	python ~/TIFseq2/combine_ends.py $outdir/${name[0]}.name.bam
#	umi_tools dedup -I $outdir/${name[0]}*_sorted.bam -S $outdir/${name[0]}_unique_UMI.bam --method cluster -L $outdir/${name[0]}_unique_UMI.log --paired --output-stats $outdir/${name[0]}_unique_UMI.stats
	python ~/TIFseq2/dedup.py -I $outdir/${name[0]}*_sorted.bam -S $outdir/${name[0]}_unique_UMI.bam --method=cluster -L $outdir/${name[0]}_unique_UMI.log --paired --spliced-is-unique
done
