#!/bin/bash 
# Add RASQUAL statistics, CADD score, conservations, TFBS, and open chromatin to input file
# The output file has the following columns: 
# 1. chr
# 2. pos-1 
# 3. pos
# 4. ref
# 5. alt
# 6. rsid 
# 7. fid
# 8. GWAS posterior
# 9. eQTL posterior
# 10. CLPP score
# 11. GWAS set (binary, in=1 or not=0)
# 12. eQTL set 
# 13. allele frequency (not MAF)
# 14. pi (effect size)
# 15. number of tested SNPs
# 16. squared correlation between prior and posterior genotypes (rSNP)
# 17. p-value
# 18. adjusted p-value (BH procedure)
# 19. eQTL rank (by p-value)
# 20. MAF 
# 21. CADD raw score
# 22. CADD phred score
# 23* number of TFBS overlaps (Pouya's Roadmap motifs)
# 24* CpG (this and all features through 25 are from CADD)
# 25* priPhCons
# 26* mamPhCons
# 27* verPhCons
# 28* priPhyloP
# 29* mamPhyloP
# 30* verPhyloP
# 31* GerpN
# 32* GerpS
# 33* GerpRS
# 34* GerpRSpval
# 35* fitCons
# 36* TFBS
# 37* TFBSPeaks
# 38* TFBSPeaksMax
# 39* ATACseq peaks for sample 2305
# 40-44* 1KG allele frequency from 5 super populations (5 columns: EAS, AMR, AFR, EUR, SAS)


in_file=$1
out_file=$2
# in_file=../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.bed
# out_file=../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.allFeatures.bed


# Add RASQUAL statistics to putative variants selected by eCAVIAR: 
bash mpra/add_RASQUAL.sh ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.bed ../processed_data/rasqual/output/ /srv/scratch/bliu2/add_RASQUAL/ ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.withRASQUAL.bed


# Add conservation, TFBS, open chromatin, CADD score to putative variants selected by eCAVIAR: 
bedtools intersect -wa -wb -loj -a ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.withRASQUAL.bed -b ../data/features/chr*.allFeatures.bed.gz > ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.withRASQUAL.withAnno.bed


# Removing redundant chrom and pos columns and reorder columns: 
cat ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.withRASQUAL.withAnno.bed | \
	awk 'BEGIN{ecaviarEnd=11; rasqaulStart=15; rasqualEnd=22; annoStart=30}{
	for(i=1; i<=ecaviarEnd; i++){ printf("%s\t",$i)}; 
	for(i=rasqaulStart; i<=rasqualEnd; i++){ printf("%s\t",$i)}; 
	for(i=annoStart; i<NF; i++){ printf("%s\t",$i)};
	print $NF}' | awk 'BEGIN{OFS="\t";print "chr","start","end","ref","alt","rsid","fid","gwasPost","eqtlPost","clpp","gwasSet","eqtlSet","af","pi","n_rsnp","rsq_rsnp","pval","padj","rank","maf","caddRaw","caddPhred","n_tfbs","CpG","priPhCons","mamPhCons","verPhCons","priPhyloP","mamPhyloP","verPhyloP","GerpN","GerpS","GerpRS","GerpRSpval","fitCons","TFBS","TFBSPeaks","TFBSPeaksMax","ATACseq","mafEAS", "mafAMR", "mafAFR", "mafEUR", "mafSAS";
	rsid_col=12;new_rsid_col=6}{
	for(i=1; i<new_rsid_col; i++){printf("%s\t",$i)};
	printf("%s\t",$rsid_col);
	for(i=new_rsid_col;i<NF;i++){if (i!=rsid_col) {printf("%s\t",$i)}}
	print $NF}'> $out_file

# rm ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.withRASQUAL.bed 
# rm ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.withRASQUAL.withAnno.bed