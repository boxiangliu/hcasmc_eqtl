# Multi-tissue eQTL mapping with metasoft
# Boxiang Liu
# 2017-12-21


## Run Metasoft (subsampled to 52, select eQTL with p<1e-3): 
# copy HCASMC eQTL data to appropriate location:
mkdir $processed_data/160816/subsampling/HCASMC
ln -s $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt.gz $processed_data/160816/subsampling/HCASMC/HCASMC_52.allpairs.txt.gz

# Concatenate eQTL result for all tissue, append tissue name, and sort according to gene ID and SNP: 
bash $scripts/160805/concatenate_eqtl_tissues.subsample.sh 

# Generate Metasoft input file:
mkdir $processed_data/160805/metasoft_input_subsample_52/ 
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py > $processed_data/160805/metasoft_input_subsample_52/metasoft_input.txt

# Split metasoft input by chromosome: 
parallel 'grep "_{}_" ../processed_data/160805/metasoft_input_subsample_52/metasoft_input.txt > ../processed_data/160805/metasoft_input_subsample_52/metasoft_input.{}.txt' ::: {1..22} X

# Run METASOFT:
mkdir $processed_data/160805/metasoft_output_subsample_52
parallel bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input_subsample_52/metasoft_input.{}.txt $processed_data/160805/metasoft_output_subsample_52/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output_subsample_52/metasoft_output.{}.mcmc.log ::: {1..22} X

# Merge Metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt

## Run Metasoft (subsample to 52, select eQTL with p<1e-2)
# generate metasoft input file:
mkdir $processed_data/160805/metasoft_input_subsample_52_p1e-2/
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py 1e-2 > $processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.txt

# split Metasoft input by chromosome: 
parallel 'grep "_{}_" ../processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.txt > ../processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.{}.txt' ::: {1..22} X

# run METASOFT:
mkdir $processed_data/160805/metasoft_output_subsample_52_p1e-2
parallel bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.{}.txt $processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.{}.mcmc.log ::: {1..22} X

# Merge Metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt

# plot heatmap of m-values: 
Rscript $scripts/160805/plot_mvalue_heatmap.R 

