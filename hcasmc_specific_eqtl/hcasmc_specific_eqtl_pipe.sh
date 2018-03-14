# setup:
mkdir hcasmc_specific_eqtl ../processed_data/hcasmc_specific_eqtl ../figures/hcasmc_specific_eqtl


# find HCASMC specific eQTLs (full sample): 
parallel -j5 Rscript ../scripts/hcasmc_specific_eqtl/find_tissue_specific_eqtls.R ../processed_data/160805/metasoft_output/metasoft_output.{2}.mcmc.txt ../figures/hcasmc_specific_eqtl/fullsample/{1}/ ../processed_data/hcasmc_specific_eqtl/fullsample/{1}/stat.{2}.txt ../processed_data/hcasmc_specific_eqtl/fullsample/{1}/tissue_specific_eqtl.{2}.txt {1} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt ::: {1..22}


# find tisue specific eQTLs (subsampled to 52):
parallel -j5 Rscript ../scripts/hcasmc_specific_eqtl/find_tissue_specific_eqtls.R ../processed_data/160805/metasoft_input_subsample_52/metasoft_output.{2}.mcmc.txt ../figures/hcasmc_specific_eqtl/subsample_52/{1}/ ../processed_data/hcasmc_specific_eqtl/subsample_52/{1}/stat.{2}.txt ../processed_data/hcasmc_specific_eqtl/subsample_52/{1}/tissue_specific_eqtl.{2}.txt {1} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt ::: {1..22}


# Plot examples of HCASMC-specific eQTL:
bash hcasmc_specific_eqtl/plot_examples_of_hcasmc_eqtl.sh


# make manhattan plot: 
Rscript hcasmc_specific_eqtl/manhattan.R


# Calculate eQTL specificity using information theory (fill NAs with 1e-8):
cat ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python $scripts/hcasmc_specific_eqtl/calc_eqtl_specificity.py ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/eqtl.1e-8.txt.gz ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/specificity.1e-8.txt.gz 1e-8


# Calculate eQTL specificity using information theory (fill NAs using mean imputation):
cat ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python $scripts/hcasmc_specific_eqtl/calc_eqtl_specificity.py ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/eqtl.mean.txt.gz ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/specificity.mean.txt.gz -1


# select HCASMC specific eQTL based on QSI, and plot distribution of their p-values: 
bash hcasmc_specific_eqtl/plot_hcasmc_specific_eqtls_pval.sh


# naively select HCASMC specific eQTL: 
cat ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python hcasmc_specific_eqtl/naive_select.py > ../processed_data/hcasmc_specific_eqtl/hcasmc_specific_eqtl.naive_select.txt


# Make forest plot: 
bash hcasmc_specific_eqtl/plot_hcasmc_specific_eqtls.sh 

