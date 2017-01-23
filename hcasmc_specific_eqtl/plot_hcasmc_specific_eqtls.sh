mkdir -p /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl

# IER3:
grep "ENSG00000137331.11_6_30507577_A_AC_b37" ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt > /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000137331.11_6_30507577_A_AC_b37.txt
Rscript hcasmc_specific_eqtl/forest_plot.R /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000137331.11_6_30507577_A_AC_b37.txt ../figures/hcasmc_specific_eqtl/ENSG00000137331.11_6_30507577_A_AC_b37.pdf
# ACTA2: 
grep "ENSG00000107796.8_10_90666952_C_T_b37" ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt > /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000107796.8_10_90666952_C_T_b37.txt
Rscript hcasmc_specific_eqtl/forest_plot.R /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000107796.8_10_90666952_C_T_b37.txt ../figures/hcasmc_specific_eqtl/ENSG00000107796.8_10_90666952_C_T_b37.pdf
# KANSL1: 
grep "ENSG00000120071.8_17_43654468_C_T_b37" ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt > /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000120071.8_17_43654468_C_T_b37.txt
Rscript hcasmc_specific_eqtl/forest_plot.R /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000120071.8_17_43654468_C_T_b37.txt ../figures/hcasmc_specific_eqtl/ENSG00000120071.8_17_43654468_C_T_b37.pdf
# MT1B: 
grep "ENSG00000169688.10_16_56666848_G_A_b37" ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt > /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000169688.10_16_56666848_G_A_b37.txt
Rscript hcasmc_specific_eqtl/forest_plot.R /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000169688.10_16_56666848_G_A_b37.txt ../figures/hcasmc_specific_eqtl/ENSG00000169688.10_16_56666848_G_A_b37.pdf
# NUPR2: 
grep "ENSG00000185290.3_7_55803970_G_A_b37" ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt > /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000185290.3_7_55803970_G_A_b37.txt
Rscript hcasmc_specific_eqtl/forest_plot.R /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl/ENSG00000185290.3_7_55803970_G_A_b37.txt ../figures/hcasmc_specific_eqtl/ENSG00000185290.3_7_55803970_G_A_b37.pdf

rm -r /srv/scratch/bliu2/HCASMC_eQTL/hcasmc_specific_eqtl