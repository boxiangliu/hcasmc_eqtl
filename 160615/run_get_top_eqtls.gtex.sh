num=$1
outdir=$2
inlist=($(ls /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6_fastQTL_FOR_QC_ONLY/*.egenes.txt.gz))
echo ${#inlist[@]} samples
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160615/
for input in ${inlist[@]}; do
	Rscript $scripts/get_top_eqtls.gtex.R \
	-input=$input \
	-outdir=$outdir \
	-num=$num
done 
