out_dir='../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/'
hapmap_fn='/srv/persistent/bliu2/shared/ldscore/w_hm3.snplist'
mkdir -p $out_dir

~/tools/ldsc/munge_sumstats.py \
	--sumstats ../data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt \
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


~/tools/ldsc/munge_sumstats.py \
	--sumstats ../data/gwas/scz/scz2.snp.results.txt \
	--out $out_dir/scz \
	--merge-alleles $hapmap_fn \
	--N-cas 36989 \
	--N-con 113075  \
	--info info \
	--info-min 0.4 \
	--p p \
	--a1 a1 \
	--a2 a2 \
	--a1-inc \
	--snp snpid \
	--signed-sumstats or,1 