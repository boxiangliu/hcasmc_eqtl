# Create a lsit of EUR samples: 
eur_sample_fn='../processed_data/gwas_atacseq_overlap/prepare_vcf/eur_sample.txt'
system(sprintf('grep EUR ../../shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > %s',eur_sample_fn))


# Prepare VCF files: 
for (c in 1:22){
	# Subset to EUR samples, change rs ID to chr:pos, and remove variants with duplicate genomic positions:
	print(sprintf('INFO - chr%s',c))
	print('INFO - preparing VCF...')
	vcf_out_fn=sprintf('%s/chr%s.vcf.gz',tmp_dir,c)
	command=sprintf("../../tools/bcftools/bcftools view -S %s ../../shared/1000genomes/phase3v5a/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Ou | ../../tools/bcftools/bcftools annotate --set-id '%%CHROM:%%POS' -Ou | ../../tools/bcftools/bcftools norm -d 'any' -Oz > %s",eur_sample_fn,c,vcf_out_fn)
	system(command)
}
