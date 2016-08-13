wd=/srv/persistent/bliu2/tools/ForestPMPlot/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures
input=$processed_data/160805/metasoft_input.txt
output=$processed_data/160805/metasoft_output.mcmc.txt
study_name=$processed_data/160805/Metasoft_tissue_order.txt
study_order=$processed_data/160805/Metasoft_tissue_idx.txt
gene_name=ENSG00000227232.4
rsid=ENSG00000227232.4_1_734349_T_C_b37
out_file=$figures/160805/$rsid.pdf 

cd $wd
python pmplot.py $input $output $study_name $study_order $rsid $gene_name $out_file
