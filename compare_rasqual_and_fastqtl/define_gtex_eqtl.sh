in_dir=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_allpairs_FOR_QC_ONLY
out_dir=../processed_data/compare_rasqual_and_fastqtl/gtex_eqtl/
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

# Append tissue column: 
echo INFO - appending tissue name...
for f in $(ls $in_dir/*allpairs.txt.gz); do
	f=$(basename $f)
	tissue=${f/_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz/}
	echo INFO - $tissue: $f
	zcat $in_dir/$f | awk -v tissue=$tissue 'BEGIN{OFS="\t"}{print $0,tissue}' > $out_dir/${f/allpairs.txt.gz/allpairs.txt}
done


# Concatenate: 
echo INFO - concatenating...
cat $out_dir/*.allpairs.txt > $out_dir/all_tissues.txt


# Sort by fid and sid:
echo INFO - sorting...
sort -k1,1 -k2,2 $out_dir/all_tissues.txt > $out_dir/all_tissues.sorted.txt


# Select eQTL with lowest p-value:
echo INFO - getting lowest p-value...
cat $out_dir/all_tissues.sorted.txt | python compare_rasqual_and_fastqtl/get_lowest_pval.py | gzip > $out_dir/all_tissues.lowest_pval.txt.gz


# Multiple hypothesis correction:
echo INFO - running multiple hypothesis correction with TreeQTL...
bash compare_rasqual_and_fastqtl/adjust_pvalue.sh \
	$out_dir/all_tissues.lowest_pval.txt.gz \
	$out_dir/tmp/ \
	$out_dir/TreeQTL/


# Remove intermediate files:
rm $out_dir/*allpairs.txt
rm $out_dir/all_tissues.txt
gzip $out_dir/all_tissues.sorted.txt

