in_dir=../processed_data/rasqual/output_merged/
tmp_dir=../processed_data/rasqual/output_merged/tmp
out_dir=$in_dir/treeQTL

if [[ ! -d $tmp_dir ]]; then mkdir -p $tmp_dir; fi
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

cat $in_dir/expressed_genes.pval.txt | grep -v "fid" | \
	awk 'BEGIN{OFS="\t"; print "SNP","gene","beta","t-stat","p-value","FDR"}{print $3"_"$4"_"$5"_"$6"_b37",$28,-1,-1,$26,-1}' | \
	sed 's/chr//' > $tmp_dir/meqtl.txt

Rscript rasqual/adjust_pvalue.R $tmp_dir/meqtl.txt $out_dir

rm -r $tmp_dir/
