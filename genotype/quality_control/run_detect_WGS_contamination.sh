#!/bin/bash
# run detect_WGS_contamination.sh
# bosh
# durga

# dir1=../processed_data/genotype/quality_control/detect_WGS_contamination/
# vim $dir1/Caucasian.txt # changed 1508 to 2109, 2999 to 289727, 317155 to 313605
# vim $dir1/Hispanic.txt # changed 1401 to CA1401, 2105 to 2102, added 1848 and 1858
# vim $dir1/AA.txt # added 24635


# paths:
scripts=genotype/quality_control/


# run for caucasian samples: 
n=0
while read sample; do
n=$((n+1))
if [[ n -gt 20 ]];then
	wait
	n=0
fi 
echo $sample
bash $scripts/detect_WGS_contamination.sh $sample EUR &
done < ../processed_data/genotype/quality_control/detect_WGS_contamination/Caucasian.txt


# Hispanic samples:
while read sample; do
echo $sample
bash $scripts/detect_WGS_contamination.sh $sample AMR &
done < ../processed_data/genotype/quality_control/detect_WGS_contamination/Hispanic.txt


# African samples:
while read sample; do
echo $sample
bash $scripts/detect_WGS_contamination.sh $sample AFR &
done < ../processed_data/genotype/quality_control/detect_WGS_contamination/AA.txt


# Asian samples:
while read sample; do
echo $sample
bash $scripts/detect_WGS_contamination.sh $sample EAS &
done < ../processed_data/genotype/quality_control/detect_WGS_contamination/Asian.txt
