out_dir=../processed_data/atacseq/qc/samtools_stats/
mkdir -p $out_dir

samtools_stats(){
	fn=$1
	basename=$(basename $fn)
	out_dir=$2
	samtools stats $fn > $out_dir/$basename.stats
}

export -f samtools_stats

parallel samtools_stats {} $out_dir ::: $(ls /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/*/out/align/*/*PE2SE.bam)