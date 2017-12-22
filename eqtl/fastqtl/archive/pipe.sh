# nominal pass to map eQTLs: 
bash $scripts/160805/run_fastqtl.wrap.sh \
	$data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	$processed_data/160530/combined.filter.norm.bed.gz \
	$processed_data/160805/covariates.pc4.peer8.gender_letter.tsv.gz \
	$processed_data/160805/hcasmc.eqtl.pc4.peer8.txt \
	"--window 1e6"

# run p-value correction:
Rscript $scripts/160629/fastqtl_nominal_pvalue_corrections.R $processed_data/160805/hcasmc.eqtl.pc4.peer8.txt $processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt
