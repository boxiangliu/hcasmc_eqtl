#!/usr/bin/env python
# bosh liu
# 2016/04/18
# durga
# re-run failed imputation jobs 

from utils import *
REFERENCE_DIR='/srv/persistent/bliu2/shared/haplotype_reference/1000G_phase3/1000GP_Phase3/'

def impute2(params, Ne = 20000, prephased = True, dry_run = False):
	umimputed = params[0]
	lb = params[1]
	ub = params[2]
	imputed_prefix = params[3]
	reference_hap = params[4]
	reference_legend = params[5]
	genetic_map = params[6]

	if prephased:
		cmd = 'impute2 -k_hap 1000 -use_prephased_g -known_haps_g %s -h %s -l %s -m %s -int %s %s -Ne %i -o %s'%(umimputed,reference_hap,reference_legend,genetic_map,lb,ub,Ne,imputed_prefix)
		call(cmd, dry_run=dry_run)
	else: 
		pass 
		# future development
	print('[ impute2 ] finished.')

wd = sys.argv[1]

setwd(wd)

cmd = 'grep -L "%Concordance" *_summary'
failed_jobs = check_output(cmd)
failed_jobs = failed_jobs.strip().split('\n')
params = list()
for line in failed_jobs:
	print line
	split_line = line.strip().split('.')
	chrom = split_line[9]
	region = split_line[12]
	split_region = region.split('_')
	lb = split_region[0]
	ub = split_region[1]
	unimputed = "recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.%s.phased.haps"%chrom
	imputed_prefix = "recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.%s.phased.imputed.%s_%s"%(chrom,lb,ub)
	reference_hap = "%s/1000GP_Phase3_%s.hap.gz"%(REFERENCE_DIR,chrom)
	reference_legend = "%s/1000GP_Phase3_%s.legend.gz"%(REFERENCE_DIR,chrom)
	genetic_map = "%s/genetic_map_%s_combined_b37.txt"%(REFERENCE_DIR,chrom)
	params.append((unimputed, lb, ub, imputed_prefix, reference_hap, reference_legend, genetic_map))

print len(params), "jobs"
pool = Pool(processes=1)
pool.map(impute2, params)
