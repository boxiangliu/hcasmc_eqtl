# Subset to top 10000 genes:
Rscript eqtl/peer/subset_top_genes.R ../processed_data/rnaseq/preprocess/combine_rpkm/combined.filter.rpkm ../processed_data/eqtl/peer/combined.filter.top10000.rpkm 10000

# quantile normalize rpkm against other samples and 
# quantile normalize rpkm for each gene: 
Rscript eqtl/peer/normalize_rpkm.R ../processed_data/eqtl/peer/combined.filter.top10000.rpkm ../processed_data/eqtl/peer/combined.filter.top10000.norm.rpkm

# transpose rpkm: 
Rscript eqtl/peer/transpose_rpkm.R ../processed_data/eqtl/peer/combined.filter.top10000.norm.rpkm ../processed_data/eqtl/peer/combined.filter.top10000.norm.t.rpkm


# The data matrix is assumed to have N rows and G columns, where N is the number of samples, and G is the number of genes:
# an example input file is at /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160527/Whole_Blood.rpkm.log2.ztrans.txt
bash eqtl/peer/GetPeerExtended.sh ../processed_data/eqtl/peer/combined.filter.top10000.norm.t.rpkm ../processed_data/eqtl/peer/


