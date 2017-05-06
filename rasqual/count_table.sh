# Prepare RASQUAL read count table: 
cat ../data/rnaseq2/read_count/rnaseqc/rnaseqc.hcasmc_eqtl.reads.gct | grep -v -e "#1.2" -e "56238" -e "Name" | cut -f1,3- > ../processed_data/rasqual/expression/Y.txt
