#!/bin/bash
dst=$1

tissues=$(ls /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_FOR_QC_ONLY/*.rpkm.gct)
for tissue in ${tissues[@]}; do
	echo "working on" $tissue
	bname=$(basename $tissue)
	cat $tissue | tail -n +3 | cut -f1,3- > $dst/${bname/.gct/}
done 
