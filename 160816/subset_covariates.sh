# subset to 10 PEER factors for eQTL calling on subsampled GTEx tissues.
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/
cov_files=($(ls /srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_All_Tissues_covariates_FOR_QC_ONLY/*_Analysis.covariates.txt))
mkdir $processed_data/160816/covariates/
for input in ${cov_files[@]}; do
	base=$(basename $input)
	grep -P "(ID|C[1-3]|InferredCov([1-9]{1}|10)|gender|Platform)\t" $input > $processed_data/160816/covariates/${base/txt/subsample.52.txt}
done 
