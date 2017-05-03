# process ATACseq data with Kundaje pipeline: 
bds ~/atac_dnase_pipelines/atac.bds -species hg19 -nth 12 -title 2305 -out_dir ~/atacseq/2305/SF/out \
 -fastq1_1 ~/atacseq/2305/SF/fastq/CA2305-SF1_S2_concat_R1_001.fastq.gz -fastq1_2  ~/atacseq/2305/SF/fastq/CA2305-SF1_S2_concat_R2_001.fastq.gz \
 -fastq2_1 ~/atacseq/2305/SF/fastq/CA2305-SF_S1_concat_R1_001.fastq.gz -fastq2_2 ~/atacseq/2305/SF/fastq/CA2305-SF_S1_concat_R2_001.fastq.gz
