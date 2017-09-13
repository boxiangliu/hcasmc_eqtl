export out_dir=../processed_data/rna_wgs_match/hg19toGRCh37/
mkdir -p $out_dir


hg19toGRCh37(){
	sample=${1/\//}
	echo $sample
	samtools view -h ../data/rnaseq2/alignments/$sample/Aligned.out.sorted.bam chr22 | \
	sed -e '/^@SQ/s/SN\:chr/SN\:/' -e '/^[^@]/s/\tchr/\t/g' | samtools view -b > $out_dir/$sample.chr22.bam
	samtools index $out_dir/$sample.chr22.bam
}

export -f hg19toGRCh37

parallel hg19toGRCh37 {} :::: ../data/rnaseq2/alignments/sample_list.txt