export sumstats_fn='../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/cad.sumstats.gz'
export baseline_annotation_dir='/srv/persistent/bliu2/shared/ldscore/baseline/'
export tissue_specific_annotation_dir='../processed_data/gwas_gene_overlap/ldscore_regression/ldscore_merged/'
export weight_dir='/srv/persistent/bliu2/shared/ldscore/weights_hm3_no_hla/'
export frq_dir='/srv/persistent/bliu2/shared/ldscore/1000G_frq/'
export out_dir='../processed_data/gwas_gene_overlap/ldscore_regression/partition_heritability_merged/'

mkdir -p $out_dir

python ~/tools/ldsc/ldsc.py \
--h2 $sumstats_fn \
--ref-ld-chr $tissue_specific_annotation_dir/merged. \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $out_dir/cad.merged \
--print-coefficients
