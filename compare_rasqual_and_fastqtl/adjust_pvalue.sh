in_fn=$1
tmp_dir=$2
out_dir=$3

if [[ ! -d $tmp_dir ]]; then mkdir -p $tmp_dir; fi
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

zcat $in_fn | \
	awk 'BEGIN{OFS="\t"; print "SNP","gene","beta","t-stat","p-value","FDR"}{print $2,$1,-1,-1,$4,-1}' \
	> $tmp_dir/meqtl.txt

Rscript rasqual/adjust_pvalue.R $tmp_dir/meqtl.txt $out_dir

rm -r $tmp_dir/
