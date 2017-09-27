export sumstats_fn='../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/cad.sumstats.gz'
export baseline_annotation_dir='/srv/persistent/bliu2/shared/ldscore/baseline/'
export tissue_specific_annotation_dir='../processed_data/gwas_gene_overlap/ldscore_regression/ldscore/'
export weight_dir='/srv/persistent/bliu2/shared/ldscore/weights_hm3_no_hla/'
export frq_dir='/srv/persistent/bliu2/shared/ldscore/1000G_frq/'
export out_dir='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability/'

mkdir -p $out_dir


# Functions:
partition_heritability(){
tissue=$1
echo $tissue

if [[ -e $out_dir/cad."$tissue".log ]]; then
echo "INFO - $tissue already processed. skipping..."
return 0
fi


python ~/tools/ldsc/ldsc.py \
--h2 $sumstats_fn \
--ref-ld-chr $tissue_specific_annotation_dir/"$tissue".,$baseline_annotation_dir/baseline. \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $out_dir/cad."$tissue" \
--print-coefficients

}

export -f partition_heritability

partition_heritability_no_baseline(){
tissue=$1
echo $tissue

if [[ -e $out_dir/cad."$tissue".nobaseline.log ]]; then
echo "INFO - $tissue already processed. skipping..."
return 0
fi


python ~/tools/ldsc/ldsc.py \
--h2 $sumstats_fn \
--ref-ld-chr $tissue_specific_annotation_dir/"$tissue". \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $out_dir/cad."$tissue".nobaseline \
--print-coefficients

}

export -f partition_heritability_no_baseline

# Get tissue names:
ls $tissue_specific_annotation_dir/*22.l2.ldscore.gz | sed "s/.22.l2.ldscore.gz//" | sed "s:$tissue_specific_annotation_dir/::" | uniq | sort > $out_dir/tissue_list.txt


# Partition heritability:
parallel -j10 partition_heritability {} :::: $out_dir/tissue_list.txt
parallel -j10 partition_heritability_no_baseline {} :::: $out_dir/tissue_list.txt
