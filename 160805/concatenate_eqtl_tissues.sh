# get all inputs: 
inputs=($(ls /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/*_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz))

# unzip each tissue and append tissue name: 
n=0
for input in ${inputs[*]}; do
	n=$(($n+1))

	# get tissue name:
	tissue=$(echo $input | grep -oP "(?<=v6p_fastQTL_allpairs_FOR_QC_ONLY/)(.+?)(?=_Analysis)")
	echo $tissue


	# unzip:
	zcat $input | awk -v tissue=$tissue 'BEGIN{OFS="\t"}{print $1"_"$2,$4,$5,$6,tissue}' > ${input/txt.gz/txt} &


	# pause after launching 10 jobs:
	if [[ $n -ge 10 ]]; then
		wait 
		n=0
	fi
done
wait


# concatenate all tissues: 
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/*_Analysis.v6p.FOR_QC_ONLY.allpairs.txt > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.txt



# sort: 
sort -S 5G --parallel 10 -T /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/tmp -k1,1 -k5,5 -V /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.txt > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt
