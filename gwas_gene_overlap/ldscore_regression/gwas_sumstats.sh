in_fn='../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt'
out_dir='../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/'
hapmap_fn='/srv/persistent/bliu2/shared/ldscore/w_hm3.snplist'
mkdir -p $out_dir

~/tools/ldsc/munge_sumstats.py \
	--sumstats $in_fn \
	--out $out_dir/cad \
	--merge-alleles $hapmap_fn \
	--N-cas 60801 \
	--N-con 123504 \
	--frq effect_allele_freq \
	--info median_info \
	--info-min 0.4 \
	--nstudy n_studies \
	--nstudy-min 29 \
	--p p_dgc \
	--a1 effect_allele \
	--a2 noneffect_allele \
	--a1-inc \
	--snp markername \
	--signed-sumstats beta,0 