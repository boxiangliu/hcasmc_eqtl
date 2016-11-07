# The Y.txt file is tab delimited and has the following format: 
# gene_name sample1_count sample2_count etc
# what the heck is K? how to generate it?
# the only thing needed is the Y.txt

# the covariates file X can be generated using provided script


# how to generate the AS field in vcf? 
# how to generate the AP field in VCF? 
# does the imputed file has rsq?
# AP, GL and GP, DS, rsq
# recalibrated_biallelic_variants.beagle.vcf has 


# joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.GRCh37.vcf.gz
# AR2=0.86;DR2=0.89;AF=0.36

# 160604_phasing/conform_gt/mod.chr10.vcf.gz does not have any INFO flag

# processed_data/160604_phasing/phased_and_imputed has AR2, DR2, AF INFO fields 

# processed_data/160604_phasing/phased/ does not have any INFO field. However, I think 
# we can assume that AR2=1 and DR2=1

# how to prepare GC content file.

# why are only the top variant get  
# generate read count table (Y file) and offset table (K.txt)


# generate VCF file with AS counts: 
# which bam files should I use? 
# Use the WASP filtered reads 

# post imputation filtering: 
# HWE P-value < 1E-6, monomorphic variants, indels of length >51 bp, and imputation quality score INFO<0.4
# 1. reheader
# 2. calculate HWE using vcftools 
# 3. filter based on HWE, imputation quality, monomorphic variants,
# 4. filter indel length using awk


# select_genes.R 
# download genes in GWAS catelog 
# 2015 nikpay paper. 
# 