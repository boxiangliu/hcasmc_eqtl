tmp_dir=../processed_data/eqtl_and_atacseq/tmp
[[ ! -d $tmp_dir ]] && mkdir -p $tmp_dir
sort -r -k2 ../processed_data/eqtl_and_atacseq/specificity.mean.txt > $tmp_dir/specificity.mean.sorted.txt
head -n10000 $tmp_dir/specificity.mean.sorted.txt | cut -f1 > $tmp_dir/specificity.mean.sorted.10000.txt
grep -f $tmp_dir/specificity.mean.sorted.10000.txt $processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt > $tmp_dir/specificity.mean.sorted.10000.eql.txt
Rscript hcasmc_specific_eqtl/plot_hcasmc_specific_eqtls_pval.R
rm -r $tmp_dir