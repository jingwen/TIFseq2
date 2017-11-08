#!/bin/bash

module load bcl2fastq
module load FastQC
module load samtools

OPTIND=1

while getopts "hd:i:o:x:" opt; do
	case $opt in
	h) 
		echo "Usage: demultiplex.sh -i <input directory> -d <sequence settings> -x 1,1 -o <output directory>"
		exit 0
		;;
	d) index=$OPTARG;;
	i) indir=$OPTARG;;
	o) outdir=$OPTARG;;
	x) idx_mis=$OPTARG;;
	esac
done

indir=${indir%/}
outdir=${outdir%/}
sampleSheet=$outdir/sampleSheet.csv
fastq=$outdir/fastq

if [ ! -d "$indir" ];then
	echo "Input directory does not exist!"
	exit 0
fi

if [ ! -d "$outdir" ];then
	echo "Output directory does not exist! Generate directory."
	mkdir $outdir
fi

./prep_sampleSheet.sh $index > $sampleSheet
bcl2fastq -R $indir -o $fastq --sample-sheet $sampleSheet --no-lane-splitting --barcode-mismatches $idx_mis 
./demultiplex_stats.sh $fastq/Stats/DemultiplexingStats.xml > $fastq/demultiplex_stat.txt
mkdir $outdir/fastqc
for i in $fastq/fastq/*001.fastq.gz; do fastqc -t -o $outdir/fastqc $i; done
mkdir $outdir/mapping
./cut_align.sh -i $fastq/fastq -o $outdir/mapping -r $ref

