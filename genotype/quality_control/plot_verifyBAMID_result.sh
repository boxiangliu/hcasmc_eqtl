# plot verifyBamID result:
dir1=../processed_data/genotype/quality_control/detect_WGS_contamination/
cat $dir1/verifyBAMID.*.selfSM | awk 'BEGIN{OFS="\t"} {if ($1!="#SEQ_ID") print $1,$7}' > $dir1/verifyBAMID.combined.tsv
Rscripts genotype/quality_control/plot_verifyBAMID_result.R