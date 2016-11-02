/srv/persistent/bliu2/tools/MetaXcan/software/MetaXcan.py \
--beta_folder ../processed_data/metaXcan/intermediate/beta \
--weight_db_path /srv/persistent/bliu2/tools/MetaXcan/software/data/GTEx-V6p-1KG-2016-08-18/TW_Artery_Coronary_0.5.db \
--covariance /srv/persistent/bliu2/tools/MetaXcan/software/data/GTEx-V6p-1KG-2016-08-18/TW_Artery_Coronary.txt.gz \
--gwas_folder /srv/persistent/bliu2/HCASMC_eQTL/data/gwas/ \
--gwas_file_pattern "cad.add.160614.website.plink.txt.gz" \
--compressed \
--beta_column BETA \
--pvalue_column P \
--output_file ../processed_data/metaXcan/results/Artery_Coronary.csv