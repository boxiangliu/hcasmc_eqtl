# get all inputs: 
inputs=($(ls /srv/persistent/bliu2/HCASMC_eQTL/processed_data//160816/subsampling/*/*_52.allpairs.txt.gz))

# unzip each tissue and append tissue name: 
n=0
for input in ${inputs[*]}; do
	n=$(($n+1))

	# get tissue name:
	tissue=$(echo $input | grep -oP "(?<=subsampling/)(.+?)(?=/)")
	echo $tissue


	# unzip:
	zcat $input | grep -v "gene_id" | awk -v tissue=$tissue 'BEGIN{OFS="\t"}{print $1"_"$2,$4,$5,$6,tissue}' > ${input/txt.gz/txt} &


	# pause after launching 10 jobs:
	if [[ $n -ge 15 ]]; then
		wait 
		n=0
	fi
done
wait


# concatenate all tissues:
echo "concatenating..."
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/*/*_52.allpairs.txt > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.txt


# sort: 
echo "sorting..."
sort -S 5G --parallel 20 -T /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/tmp -k1,1 -k5,5 -V /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.txt > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt