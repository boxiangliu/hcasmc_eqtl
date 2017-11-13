# Calculate enrichment of genomic annotation 
# using qtlBHM. 

out_dir=../processed_data/sqtl/fastQTL/enrichment/
mkdir -p $out_dir

# Functions:
# Prepare summary statistics: 
prepare_sumstat(){
	in_fn=$1
	out_fn=$2
	echo INFO - input: $in_fn
	echo INFO - output: $out_fn

	zcat $in_fn | \
	awk -F'[\t_]' 'BEGIN{OFS="\t"}{print $1"_"$2,"chr"$3"."$4,$10,$11}' | \
	bgzip > $out_fn
}

export -f prepare_sumstat

# echo INFO - preparing summary statistics
# parallel -j10 prepare_sumstat \
# ../processed_data/sqtl/fastQTL/nominal_sid/chr{}.nominal.txt.gz \
# $out_dir/chr{}.statistics.txt.gz ::: {1..22}

echo INFO - estimating enrichment
python /srv/persistent/bliu2/tools/qtlBHM/infer_causal_variants.py \
--output_prefix $out_dir/chr22_qtlBHM \
../processed_data/sqtl/fastQTL/enrichment/chr22.statistics.txt.gz \
../processed_data/sqtl/fastQTL/snpEff//all.3.bed.gz


# echo INFO - estimating enrichment
# python /srv/persistent/bliu2/tools/qtlBHM/infer_causal_variants.py \
# --output_prefix ../processed_data/sqtl/fastQTL/test_qtlBHM/test_qtlBHM \
# ../processed_data/sqtl/fastQTL/test_qtlBHM/statistics.txt.gz \
# ../processed_data/sqtl/fastQTL/test_qtlBHM/annotation.bed.gz
