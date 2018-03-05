#!/bin/bash

module load bcl2fastq
module load samtools

OPTIND=1

while getopts "hd:i:o:x:p:" opt; do
	case $opt in
	h) 
		echo "Usage: main.sh -d <sequence index info> -i <input directory> -o <output directory> -x 2,1 -p <polyA tail size>"
		exit 0
		;;
	d) index=$OPTARG;;
	i) indir=$OPTARG;;
	o) outdir=$OPTARG;;
	x) idx_mis=$OPTARG;;
	p) polyA=$OPTARG;;
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

prep_sampleSheet.awk $index > $sampleSheet
bcl2fastq -R $indir -o $fastq --sample-sheet $sampleSheet --no-lane-splitting --barcode-mismatches $idx_mis
demultiplex_stats.awk $fastq/Stats/DemultiplexingStats.xml > $fastq/demultiplex_stat.txt
preprocess.sh -I $fastq -O $outdir -A $polyA
