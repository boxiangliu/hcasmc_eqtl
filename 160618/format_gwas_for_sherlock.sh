gwas_file=$1
output_file=$2
# gwas_file=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt
# output_file=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160618/sherlock/cad.add.sherlock.txt
cat $gwas_file | awk 'BEGIN {OFS="\t"} {print $1,$11}' | grep -v "markername" > $output_file