export bam_dir=../rna_wgs_match/hg19toGRCh37/
export vcf_fn=../processed_data/subset_vcf_to_coding_region/chr22.vcf.gz
export out_dir=../processed_data/verifyBamID/
mkdir -p $out_dir

run_verifyBamID(){
	sample=${1/\//}
	echo $sample
	verifyBamID --vcf $vcf_fn \
	--bam $bam_dir/$sample.chr22.bam \
	--smID $sample \
	--out $out_dir/$sample \
	--best \
	--ignoreRG
}

export -f run_verifyBamID

parallel -j20 run_verifyBamID :::: ../data/rnaseq2/alignments/sample_list.txt