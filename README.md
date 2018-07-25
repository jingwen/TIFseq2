# TIFseq2

### Prepare sample sheet for demultiplexing
prep_sampleSheet.awk <index_info> > <sample_sheet>
### Demultiplexing
bcl2fastq -R <input_dir> -o <fastq_dir> --sample-sheet <sample_sheet> --no-lane-splitting --barcode-mismatches <mismatch>
demultiplex_stats.awk <fastq_dir>/Stats/DemultiplexingStats.xml > <fastq_dir>/demultiplex_stat.txt
### Preprocess
preprocess.sh -I <fastq_dir> -O <output_dir> -j <thread_number> -A <polyA_length>
### Alignment of 5'end and 3'end reads
STAR_align.sh -R <STAR_index_dir> -A <splicing_junction_gtf> -I <output_dir>/cutPolyA -O <STAR_output_dir> -p <thread_number> -j <max_intron_size> -m 0
