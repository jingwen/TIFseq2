# TIF-Seq2

The analysis pipeline and downstream analysis in jingwen/TIFseq2 are used for publication [TIF-Seq2 disentangles overlapping isoforms in complex human transcriptomes](https://doi.org/10.1093/nar/gkaa691). The scripts for TIF-Seq2 data pre-processing and alignment in other TIF-Seq2 related publication are in [PelechanoLab/TIFseq2](https://github.com/PelechanoLab/TIFseq2).

### Prepare sample sheet for demultiplexing
prep_sampleSheet.awk <index_info> > <sample_sheet>
### Demultiplex
bcl2fastq -R <input_dir> -o <fastq_dir> --sample-sheet <sample_sheet> --no-lane-splitting --barcode-mismatches <mismatch>
demultiplex_stats.awk <fastq_dir>/Stats/DemultiplexingStats.xml > <fastq_dir>/demultiplex_stat.txt
### Preprocess
preprocess.sh -I <fastq_dir> -O <output_dir> -j <thread_number> -A <polyA_length>
### Align TIF-Seq2 reads and remove PCR duplicates 
STAR_align.sh -R <STAR_index_dir> -A <splicing_junction_gtf> -I <output_dir>/cutPolyA -O <STAR_output_dir> -p <thread_number> -j <max_intron_size> -m 0
### Fetch boundaries of transcription isoforms (TIFs)
python boundary.py <input_bam>
### Filter internal priming of 3'end
Rscript clean_As.R <3end_ctss_path>
### Cluster 5' ends and 3' ends respectively
Rscript cluster_end.R <TIF_5end_path> <TIF_3end_path> <3'T-fill_3end_path>
### Construct transcription isoform (TIF) boundaries
Rscript form_TIF.R
