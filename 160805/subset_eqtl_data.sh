input_files=($(ls /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/*allpairs.txt.gz))
for input_file in ${input_files[@]}; do
	output_file=${input_file/txt.gz/head2000.txt}
	echo "zcat $input_file | head -n2000 > $output_file"
	zcat $input_file | head -n2000 > $output_file
done 