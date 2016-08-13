processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/
for chr in {1..22} X; do
	grep "_${chr}_" $processed_data/160805/metasoft_input/metasoft_input.txt > $processed_data/160805/metasoft_input/metasoft_input.$chr.txt
done 
