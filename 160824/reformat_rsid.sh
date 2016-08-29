n=0
for in_file in $(ls /srv/persistent/bliu2/shared/1000genomes/phase1v3/ALL.chr{1..22}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz);do
	n=$((n+1))
	bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT\_b37' -Oz -o ${in_file/vcf.gz/rsid.vcf.gz} $in_file &
	if [[ n -ge 11 ]]; then 
		wait 
		n=0
	fi
done