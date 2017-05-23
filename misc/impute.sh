#!/bin/bash 
# bosh liu
# 2016/04/10

phased=$1
reference_hap=$2
reference_legend=$3
genetic_map=$4
imputed=$5
size=$6
num_processes=0
for lb in $(seq -f %1.0f 1 5e6 $size); do
	ub=$(($lb+5000000-1))
	if [[ ub -ge $size ]]; then 
		ub=$size
	fi
	echo "lower bound:" $lb
	echo "upper bound:" $ub
	impute2 -k_hap 1000 -use_prephased_g -known_haps_g $phased -h $reference_hap -l $reference_legend -m $genetic_map -int $lb $ub -Ne 20000 -o $imputed.${lb}_${ub} &
	
	# limit the number of processes to be less than or equal 20:
	num_processes=$((num_processes+1))
	if [[ $num_processes -ge 10 ]]; then
		wait 
		num_processes=0
	fi 
done 
