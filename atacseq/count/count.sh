out_dir=../data/atacseq/fbs/count/
mkdir -p $out_dir

parallel bedtools coverage \
-a ../processed_data/atacseq/count/merged_peak.bed \
-b ../data/atacseq/fbs/wasp/rmdup/{}.keep.merge.rmdup.sort.bam '>' $out_dir/{}.count ::: 1346  1508  1522  200212  2108  2305  2356  2510  2989
