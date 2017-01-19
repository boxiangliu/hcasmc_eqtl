in_file=$1
ecaviar_dir=$2
tmp_dir=$3
out_file=$4

# in_file=../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_genes.txt
# ecaviar_dir=../processed_data/eCAVIAR/eCAVIAR_output_rasqual_111616/
# tmp_dir=/srv/scratch/bliu2/concat_eCAVIAR_variants/
# out_file=../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.txt

if [[ ! -d $tmp_dir ]]; then mkdir $tmp_dir; fi

while read line; do
	if [[ $line == "fid"* ]]; then continue ; fi
	split_line=($line)
	if [[ ${split_line[1]} > 0.01 ]]; then 
		fid=${split_line[0]}
		echo $fid
		join -j 1 -o 1.1,1.3,2.3 <(sort -k1 $ecaviar_dir/$fid.ecaviar_1_post) <(sort -k1 $ecaviar_dir/$fid.ecaviar_2_post) > $tmp_dir/$fid.post
		join -j 1 -o 1.1,1.2,1.3,2.2 <(sort -k1 $tmp_dir/$fid.post) <(sort -k1 $ecaviar_dir/$fid.ecaviar_col) > $tmp_dir/$fid.post.col
		grep -f $ecaviar_dir/$fid.ecaviar_1_set $tmp_dir/$fid.post.col | awk 'BEGIN{OFS="\t"}{print $0,"gwas"}' > $tmp_dir/$fid.post.col.gwas
		grep -f $ecaviar_dir/$fid.ecaviar_2_set $tmp_dir/$fid.post.col | awk 'BEGIN{OFS="\t"}{print $0,"eqtl"}' > $tmp_dir/$fid.post.col.eqtl
		cat $tmp_dir/$fid.post.col.gwas $tmp_dir/$fid.post.col.eqtl > $tmp_dir/$fid.long
		Rscript mpra/reshape.R $tmp_dir/$fid.long $tmp_dir/$fid.wide
		awk -v fid=$fid 'BEGIN{FS="\t|_";OFS="\t"}{print "chr"$1,$2-1,$2,$3,$4,fid,$6,$7,$8,$9,$10}' $tmp_dir/$fid.wide > $tmp_dir/$fid.wide.withGene.bed
	fi
done < $in_file
cat $tmp_dir/*.wide.withGene.bed > $out_file
rm -r $tmp_dir