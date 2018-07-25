#!/bin/bash
module load bioinfo-tools
module load FastQC
module load cutadapt

while getopts "hI:O:j:A:" opt; do
	case $opt in
	h)
		echo "Usage: preprocess.sh -I <input fastq dir> -O <output fastq dir> -j 4 -A <polyA size>"
		exit 0
		;;
	I)indir=$OPTARG;;
	O)outdir=$OPTARG;;
	j)thread=$OPTARG;;
	A)polyA=$OPTARG;;
	esac
done

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
bothdir=$outdir/bothEnds
UMIdir=$outdir/UMI
cutdir=$outdir/cutAdapter
polyAdir=$outdir/cutPolyA
if [ ! -d "$bothdir" ];then mkdir $bothdir; fi
if [ ! -d "$UMIdir" ];then mkdir $UMIdir; fi
if [ ! -d "$cutdir" ];then mkdir $cutdir; fi
if [ ! -d "$polyAdir" ];then mkdir $polyAdir; fi
if [[ $polyA < 100 ]];then error_rate=0.1; else error_rate=0.05; fi
for fq1 in $indir/*5_S*R1_001.fastq.gz; do
	fq2=${fq1/R1_001/R2_001}
	if [ ! -f $fq2 ];then
		echo "Cannot find paired read file "$fq2" ..."
		exit 0
	fi
	IFS=_
	sample=(${fq1##*/})
	unset IFS
	if [ ${sample[0]} = ${sample[2]} ];then
		if [ "${sample[1]}" -eq 3 ] && [ "${sample[3]}" -eq 5 ];then
			out3=$bothdir/${sample[0]}_3end.fastq.gz
			out5=$bothdir/${sample[0]}_5end.fastq.gz
			umi3=$UMIdir/${sample[0]}_3UMI.fastq.gz
			umi5=$UMIdir/${sample[0]}_5UMI.fastq.gz
			cut3=$cutdir/${sample[0]}_3cutadapt.fastq.gz
			cut5=$cutdir/${sample[0]}_5cutadapt.fastq.gz
			polyA3=$polyAdir/${sample[0]}_3cutPolyA.fastq.gz
			polyA5=$polyAdir/${sample[0]}_5cutPolyA.fastq.gz
			echo "Processing "${sample[0]}
			cat $fq2 $indir/"${sample[0]}"_5_"${sample[2]}"_3_S*_R1_001.fastq.gz > $out5
			cat $fq1 $indir/"${sample[0]}"_5_"${sample[2]}"_3_S*_R2_001.fastq.gz > $out3
			echo "Cuting adapters..."
			cutadapt -j $thread -a "AGGTGACCGG" -A "AGGTGACCGG" -a "AGATCGGAAG" -A "AGATCGGAAG" --nextseq-trim=20 --match-read-wildcards --minimum-length 28 -o $cut5 -p $cut3 $out5 $out3 > $cutdir/${sample[0]}_cutAdapter.log
			echo -e "Finish cutting adapters. Now extracting UMIs..."
			umi_tools extract -I $cut5 --extract-method=string --bc-pattern=NNNNNNNN --read2-in=$cut3 --stdout=$umi5 --read2-out=$umi3 --log=$UMIdir/${sample[0]}_UMI.log
#			fastqc -o $UMIdir $umi5 $umi3
			echo -e "Finish extracting UMIs. Now cutting As stretch..."
			cutadapt -j $thread -a "A{$[polyA-8]}" -G "T{$polyA}" --nextseq-trim=20 -O 1 --error-rate $error_rate --match-read-wildcards --minimum-length 20 -o $polyA5 -p $polyA3 $umi5 $umi3 > $polyAdir/${sample[0]}_cutPolyA.log
			echo "Finish processing "${sample[0]}
		fi
	fi
done
