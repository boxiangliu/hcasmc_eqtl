in_file=$1
out_file=$2
chr=$3
cat $in_file | awk -v chr=$3 '{if ($2==chr) {print $0}}' > $out_file