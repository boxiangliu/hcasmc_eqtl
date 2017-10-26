unset DISPLAY XAUTHORITY

# Variables:
export plink_dir=/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/
export hapmap_dir=/srv/persistent/bliu2/shared/ldscore/hapmap3_snps/

annot_dir=${1-"/srv/persistent/bliu2/HCASMC_eQTL/processed_data/gwas_atacseq_overlap/ldscore_regression/tissue_specific_snp_annotation/"}
out_dir=${2-"../processed_data/gwas_atacseq_overlap/ldscore_regression/ldscore_merged/"}
log_dir=${3-"../logs/gwas_atacseq_overlap/ldscore_regression/ldscore_merged/"}
echo $out_dir $log_dir
mkdir -p $out_dir $log_dir

# Functions: 
calc_ldscore(){
tissue=$1
in_dir=$2
output_dir=$3
chr=$4

echo INFO - $tissue
echo INFO - $chr

if [[ -e "$output_dir/$tissue.$chr.log" ]]; then
check_line=$(tail -n11 "$output_dir/$tissue.$chr.log" | head -n1)

if [[ $check_line == "Summary of Annotation Matrix Row Sums" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi

fi 

python ~/tools/ldsc/ldsc.py \
--l2 \
--bfile $plink_dir/1000G.mac5eur.$chr \
--ld-wind-cm 1 \
--annot "$in_dir/$tissue.$chr.annot.gz" \
--out "$output_dir/$tissue.$chr" \
--print-snps $hapmap_dir/hm.$chr.snp

}

export -f calc_ldscore

# Calculate LD score:
for dir1 in jaccard_similarity_0.3 jaccard_similarity_0.4 jaccard_similarity_0.5 all_tissue; do
for chr in {1..22}; do
echo INFO - dir1
echo INFO - chr

mkdir -p $out_dir/$dir1 $log_dir/$dir1
calc_ldscore merged $annot_dir/$dir1/ $out_dir/$dir1/ $chr
ln $annot_dir/$dir1/* $out_dir/$dir1/

done 
done


# for dir1 in all_tissue jaccard_similarity_0.5 jaccard_similarity_0.4 jaccard_similarity_0.3; do
# mkdir -p $out_dir/$dir1 $log_dir/$dir1

# parallel -j15 calc_ldscore merged $annot_dir/$dir1/ $out_dir/$dir1/ {} ::: {1..22}
# ln $annot_dir/$dir1/* $out_dir/$dir1/

# done
