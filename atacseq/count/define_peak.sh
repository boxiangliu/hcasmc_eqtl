out_dir=../processed_data/atacseq/count/
mkdir -p $out_dir

parallel zcat ../data/atacseq/fbs/{}/out/peak/macs2/overlap/*.narrowPeak.gz "|" sort -k1,1 -k2,2n "|" bgzip ">" ../data/atacseq/fbs/{}/out/peak/macs2/overlap/{}.sorted.narrowPeak.gz ::: 1346  1508  1522  200212  2108  2305  2356  2510  2989

bedtools intersect -wa -wb \
-a ../data/atacseq/fbs/2305/out/peak/macs2/overlap/*.sorted.narrowPeak.gz \
-b ../data/atacseq/fbs/{1346,1508,1522,200212,2108,2305,2356,2510,2989}/out/peak/macs2/overlap/*.sorted.narrowPeak.gz > $out_dir/intersect.bed

Rscript atacseq/count/define_peak.R