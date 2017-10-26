in_dir=../processed_data/gwas_atacseq_overlap/ldscore_regression/tissue_specific_snp_annotation/
for dir in all_tissue jaccard_similarity_0.3 jaccard_similarity_0.4 jaccard_similarity_0.5; do
for chr in {1..22}; do
echo $dir
echo $chr
gunzip $in_dir/$dir/merged.$chr.annot.gz
sed -i '1s/ /_/g' $in_dir/$dir/merged.$chr.annot
gzip $in_dir/$dir/merged.$chr.annot
done
done