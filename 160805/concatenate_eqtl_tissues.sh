# get all inputs: 
inputs=($(ls /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/*_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz))


# sort each input: 
n=0
for input in ${inputs[*]}; do
	n=$(($n+1))

	# get tissue name:
	tissue=$(echo $input | grep -oP "(?<=v6p_fastQTL_allpairs_FOR_QC_ONLY/)(.+?)(?=_Analysis)")
	echo $tissue


	# sort
	zcat $input | awk 'BEGIN{OFS="\t"}{print $1"_"$2,$4,$5,$6,"$tissue"}' | sort -k1 -V > ${input/txt.gz/sorted.txt} &


	# pause after launching 10 jobs:
	if [[ $n -ge 10 ]]; then
		wait 
		n=0
	fi
done
wait

# merge sort all input: 
sort -m /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/*_Analysis.v6p.FOR_QC_ONLY.allpairs.2.sorted.txt > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/All_Tissues_Analysis.v6p.FOR_QC_ONLY.allpairs.2.sorted.txt
