input=$1
output=$2
rm $output
while read in_file; do
	echo $in_file
	cat ../processed_data/rasqual/output/$in_file.pval.txt | \
	awk -v in_file=$in_file 'BEGIN{OFS="\t"} {if ($1!="fid") {print in_file,$3,$4,$5,$6,".",$26,$12,"."}}' | 
	sed 's/chr//' >> \
	$output
done < $input

