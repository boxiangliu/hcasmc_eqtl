#!/bin/bash
# run detect_sample_contamination_model_based.sh
# bosh
# durga

# paths:
scripts=./160526


# run for caucasian samples: 
n=0
while read sample; do
n=$((n+1))
if [[ n -gt 20 ]];then
	wait
	n=0
fi 
echo $sample
bash $scripts/detect_sample_contamination_model_based.sh $sample EUR &
done < ../processed_data/160526/detect_sample_contamination_model_based/Caucasian.txt


# Hispanic samples:
while read sample; do
echo $sample
bash $scripts/detect_sample_contamination_model_based.sh $sample AMR &
done < ../processed_data/160526/detect_sample_contamination_model_based/Hispanic.txt


# African samples:
while read sample; do
echo $sample
bash $scripts/detect_sample_contamination_model_based.sh $sample AFR &
done < ../processed_data/160526/detect_sample_contamination_model_based/AA.txt


# Asian samples:
while read sample; do
echo $sample
bash $scripts/detect_sample_contamination_model_based.sh $sample EAS &
done < ../processed_data/160526/detect_sample_contamination_model_based/Asian.txt
