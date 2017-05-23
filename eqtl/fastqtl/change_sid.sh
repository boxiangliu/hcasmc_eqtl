in_dir='../data/joint3/asvcf/'
out_dir='../data/joint3/asvcf_sid/'

if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi
for f in $(ls $in_dir/*vcf.new.gz); do
	vcf_fn=$(basename $f)
	echo $vcf_fn
	zcat $in_dir/$vcf_fn | awk 'BEGIN{OFS="\t"}{if ($1 !~ /#/) {$3=$1"_"$2"_"$4"_"$5"_b37"; gsub(/chr/,"",$3); print $0} else {print $0}}' | bgzip > $out_dir/$vcf_fn
	tabix -p vcf $out_dir/$vcf_fn
done 
