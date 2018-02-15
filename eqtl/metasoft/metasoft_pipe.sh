# Multi-tissue eQTL mapping with metasoft
# Boxiang Liu
# 2017-12-21

# Select top eSNP per gene:
Rscript eqtl/metasoft/top_eSNP.R
bash eqtl/metasoft/gen_metasoft_input.sh
bash eqtl/metasoft/metasoft.sh
bash eqtl/metasoft/find_tissue_specific_eqtls.R

