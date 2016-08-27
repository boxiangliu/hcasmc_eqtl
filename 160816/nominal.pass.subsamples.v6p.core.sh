#!/usr/bin/bash
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts
# Core bash script to run FastQTL nominal pass for a given tissue at each of its possible subsample sizes
tissue=${1}

# Run nominal pass for each subsample size
sampleSize=52


# Nominal pass of FastQTL
subsample=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/$tissue/${tissue}_subsample.52.txt
/usr/bin/python3 $scripts/160816/run_FastQTL_threaded_subsample.py /users/joed3/GTExCisEqtls/data/subsampling/VCF/gtex.v6.allchr.impute.info04.maf01.hwep1e6.constrvarids.vcf.gz \
	/users/joed3/GTExCisEqtls/data/subsampling/${tissue}/${tissue}.normalized.expr.bed.gz \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/${tissue}/${tissue}_${sampleSize} \
	--covariates /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/covariates/${tissue}_Analysis.covariates.subsample.52.txt \
	--window 1e6 --ma_sample_threshold 0 --maf_threshold 0.01 --chunks 100 --threads 10 --include_samples ${subsample}

