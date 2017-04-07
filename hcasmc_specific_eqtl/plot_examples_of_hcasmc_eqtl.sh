# make P-M plot for ACTA2: 
grep "ENSG00000107796.8_10_90666952_C_T_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.10.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/hcasmc_specific_eqtl/subsample52_examples/acta2.size52.pmplot.pdf none

# make P-M plot for IER3:
grep "ENSG00000137331.11_6_30507577_A_AC_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.6.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/hcasmc_specific_eqtl/subsample52_examples/ier3.size52.pmplot.pdf none


# make P-M plot for MT1B:
grep "ENSG00000169688.10_16_56666848_G_A_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.16.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/hcasmc_specific_eqtl/subsample52_examples/mt1b.size52.pmplot.pdf none


# make P-M plot for KANSL1:
grep "ENSG00000120071.8_17_43654468_C_T_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.17.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/hcasmc_specific_eqtl/subsample52_examples/kansl1.size52.pmplot.pdf none


# make P-M plot for NUPR2:
grep "ENSG00000185290.3_7_55803970_G_A_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.7.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/hcasmc_specific_eqtl/subsample52_examples/nupr2.size52.pmplot.pdf none
 