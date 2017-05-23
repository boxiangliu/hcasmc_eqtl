in_dir=../processed_data/rasqual/output_pval/
out_dir=../processed_data/rasqual/output_merged/
tmp_dir=../processed_data/rasqual/tmp/

if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi
if [[ ! -d $tmp_dir ]]; then mkdir -p $tmp_dir; fi

cat ../processed_data/rnaseq/preprocess/combine_rpkm/combined.filter.rpkm | cut -f1 | grep -v Name > $out_dir/expressed_genes.txt

head -n1 ../processed_data/rasqual/output_pval/chr1/ENSG00000000457.9_SCYL3.pval.txt | awk 'BEGIN{OFS="\t"}{print $0, "gene_id"}' > $out_dir/expressed_genes.pval.txt


for i in `seq 22`; do
	echo INFO - linking chr$i
	ln $in_dir/chr$i/* $tmp_dir
done



n=0
while read gene_id; do
	n=$((n+1))
	echo INFO - $n: $gene_id
	cat $tmp_dir/$gene_id*pval.txt | grep -v 'fid' | awk -v gene_id=$gene_id 'BEGIN{OFS="\t"}{print $0, gene_id}' >> $out_dir/expressed_genes.pval.txt
done < $out_dir/expressed_genes.txt

rm -r $tmp_dir
