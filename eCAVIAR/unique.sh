in_dir=$1
out_dir=$2

# in_dir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/eCAVIAR_input2/HCASMC
# out_dir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/eCAVIAR_input3/HCASMC

if [[ ! -d $out_dir ]];then mkdir -p $out_dir;fi

# take the unique entries for id files: 
in_files=($(ls $in_dir/*.id))
for in_file in ${in_files[*]}; do
	echo "input: $in_file"
	base=$(basename $in_file)
	out_file=${base/.id/.uniq.id}
	uniq $in_file > $out_dir/$out_file
done


# gwas zscore files: 
in_files=($(ls $in_dir/*.gwas.zscore))
for in_file in ${in_files[*]}; do
	echo "input: $in_file"
	base=$(basename $in_file)
	out_file=${base/.gwas.zscore/.uniq.gwas.zscore}
	uniq $in_file > $out_dir/$out_file
done


# eqtl zscore files: 
in_files=($(ls $in_dir/*.eqtl.zscore))
for in_file in ${in_files[*]}; do
	echo "input: $in_file"
	base=$(basename $in_file)
	out_file=${base/.eqtl.zscore/.uniq.eqtl.zscore}
	uniq $in_file > $out_dir/$out_file
done
