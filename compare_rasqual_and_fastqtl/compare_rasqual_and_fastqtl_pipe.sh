# Setup: 
mkdir compare_rasqual_and_fastqtl ../processed_data/compare_rasqual_and_fastqtl

# Subset the gene, SNP, p-value, q-value, and effect size column of fastQTL and RASQUAL result: 
bash compare_rasqual_and_fastqtl/subset.sh

# Calculate the percentage of RASQUAL eQTL also discovered fastQTL: 
bash compare_rasqual_and_fastqtl/compare.R

# Replication in GTEx:
bash compare_rasqual_and_fastqtl/define_gtex_eqtl.sh
bash compare_rasqual_and_fastqtl/gtex_replication.R

# Overlap with ATACseq:
bash compare_rasqual_and_fastqtl/top_snp_per_gene.R
bash compare_rasqual_and_fastqtl/atacseq_overlap.R

# Compare eGenes:
bash compare_rasqual_and_fastqtl/compare_eGenes.R