# eQTL analysis pipeline
# 2017-11-22
# Boxiang Liu

# Setup:
mkdir eqtl ../processed_data/eqtl ../figures/eqtl

# extract PEER factors:
bash eqtl/peer/peer_factor.sh
Rscript eqtl/peer/peer_factor_correlation.R


# Matrix eQTL:
bash eqtl/matrix_eqtl/combine_covariates.sh
Rscript eqtl/matrix_eqtl/covariates_correlation.R
bash eqtl/matrix_eqtl/matrix_eqtl.sh


# FastQTL:
bash eqtl/fastqtl/change_sid.sh
bash eqtl/fastqtl/fastqtl.nominal.sh
bash eqtl/fastqtl/adjust_pvalue.nominal.sh


# Enrichment analysis:
bash eqtl/enrichment/vsea.R