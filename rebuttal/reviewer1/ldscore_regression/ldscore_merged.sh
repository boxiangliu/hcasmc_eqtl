unset DISPLAY XAUTHORITY

# Variables:
export annot_dir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rebuttal/reviewer1/ldscore_regression/tissue_specific_snp_annotation/
export plink_dir=/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/
export hapmap_dir=/srv/persistent/bliu2/shared/ldscore/hapmap3_snps/
export out_dir=../processed_data/rebuttal/reviewer1/ldscore_regression/ldscore_merged/
export log_dir=../logs/rebuttal/reviewer1/ldscore_regression/ldscore_merged/
mkdir -p $out_dir $log_dir

# Functions: 
calc_ldscore(){
tissue=$1
annot_dir=$2
out_dir=$3
chr=$4

echo INFO - $tissue
echo INFO - $chr

if [[ -e "$out_dir/$tissue.$chr.log" ]]; then
check_line=$(tail -n11 "$out_dir/$tissue.$chr.log" | head -n1)

if [[ $check_line == "Summary of Annotation Matrix Row Sums" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi

fi 

python ~/tools/ldsc/ldsc.py \
--l2 \
--bfile $plink_dir/1000G.mac5eur.$chr \
--ld-wind-cm 1 \
--annot "$annot_dir/$tissue.$chr.annot.gz" \
--out "$out_dir/$tissue.$chr" \
--print-snps $hapmap_dir/hm.$chr.snp

}

export -f calc_ldscore

# Calculate LD score:
# for dir1 in tissue_specific_gene tissue_specific_gene_no_sm tissue_specific_gene_no_sm_no_blood; do
for dir1 in tissue_specific_gene_no_sm tissue_specific_gene_no_sm_no_blood; do
for dir2 in top200 top500 top1000 top2000 top4000 top8000; do
mkdir -p $out_dir/$dir1/$dir2
parallel -j15 --joblog $log_dir/calc_ldscore.log calc_ldscore merged $annot_dir/$dir1/$dir2 $out_dir/$dir1/$dir2 {} ::: {1..22}
ln $annot_dir/$dir1/$dir2/* $out_dir/$dir1/$dir2
done
done

