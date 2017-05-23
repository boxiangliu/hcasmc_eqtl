# Calculate p-value:
in_dir=../processed_data/rasqual/output/
out_dir=../processed_data/rasqual/output_pval/
tmp_dir=../processed_data/rasqual/tmp

if [[ ! -d $tmp_dir ]]; then mkdir $tmp_dir; fi


for i in `seq 22`; do
	echo INFO - chr$i
	ls $in_dir/chr$i | sed 's/.txt//' > $tmp_dir/chr$i.tmp
	if [[ ! -d $out_dir/chr$i ]]; then mkdir -p $out_dir/chr$i; fi
	parallel -j10 Rscript rasqual/calc_pval.R $in_dir/chr$i/{}.txt $out_dir/chr$i/{}.pval.txt :::: $tmp_dir/chr$i.tmp
done 
rm -r $tmp_dir
