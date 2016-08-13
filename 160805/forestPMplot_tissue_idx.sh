processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data
tissue_idx_file=$processed_data/160805/Metasoft_tissue_idx.txt
touch $tissue_idx_file
for i in {1..45}; do 
	echo $i >> $tissue_idx_file
done 