add_RASQUAL(){
	in_file=$1
	chr=$2
	pos=$3
	grep -P "$chr\t$pos" $in_file 
}

export -f add_RASQUAL

variant_file=$1
RASQUAL_dir=$2
tmp_dir=$3
out_file=$4
if [[ -e $out_file ]]; then rm $out_file; fi
if [[ ! -d $tmp_dir ]]; then mkdir $tmp_dir; fi

while read line; do
	split_line=($line)
	fid=${split_line[5]}
	chr=${split_line[0]}
	pos=${split_line[2]}
	echo $fid $chr $pos
	in_file=$RASQUAL_dir/$fid.pval.txt
	RASQUAL_line=$(add_RASQUAL $in_file $chr $pos)
	if [[ $(echo $RASQUAL_line | wc -l) -gt 1 ]]; then 
		echo "more than 1 line found for $fid $chr $pos in $in_file"; break
	elif [[ $(echo $RASQUAL_line | wc -l) -lt 1 ]]; then
		echo "less than 1 line found for $fid $chr $pos in $in_file"; break
	else
		echo $RASQUAL_line | awk 'BEGIN{OFS="\t"}{print $3,$4-1,$4,$2,$7,$12,$18,$25,$26,exp(log(10)*$10),$27}'>> $tmp_dir/RASQUAL_stat.bed
	fi
done < $variant_file
paste $variant_file $tmp_dir/RASQUAL_stat.bed > $out_file
rm -r $tmp_dir