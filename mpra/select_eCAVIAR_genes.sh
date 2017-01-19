in_dir=$1
out_file=$2
# in_dir=../processed_data/eCAVIAR/eCAVIAR_output_rasqual_110716/
# out_file=../processed_data/mpra/ecaviar_result.txt
if [[ -e $out_file ]]; then echo "removing $out_file!"; rm $out_file; fi
echo -e "fid\tclppSum\tgwasSet\teqtlSet\tintersectSet\tunionSet" > $out_file
for col_file in $(ls $in_dir/*col); do
	base=$(basename $col_file)
	base=${base/.ecaviar_col/}
	echo $base
	clpp_sum=$(cat $col_file | cut -f2 | paste -sd+ | sed -e 's/[e]+*/\*10\^/g' | bc)
	gwas_set_card=$(wc -l $in_dir/$base.ecaviar_1_set | cut -d" " -f1)
	eqtl_set_card=$(wc -l $in_dir/$base.ecaviar_2_set | cut -d" " -f1)
	intersect_card=$(comm -12 $in_dir/$base.ecaviar_1_set $in_dir/$base.ecaviar_2_set | wc -l)
	union_card=$(comm $in_dir/$base.ecaviar_1_set $in_dir/$base.ecaviar_2_set | wc -l)
	echo -e "$base\t$clpp_sum\t$gwas_set_card\t$eqtl_set_card\t$intersect_card\t$union_card" >> $out_file
done 