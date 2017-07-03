map1_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/wasp/map1/
map1_sort_dir=/srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/wasp/map1_sort/
mkdir -p $map1_dir $map1_sort_dir
for s in 1346 1508 1522 200212 2108 2305 2356 2510 2989; do
	echo $s
	ln -s /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/$s/out/align/rep1/*PE2SE.bam $map1_dir/$s.bam
	ln -s /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/$s/out/align/rep1/*PE2SE.bam $map1_sort_dir/$s.bam
	ln -s /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/$s/out/align/rep1/*PE2SE.bam.bai $map1_sort_dir/$s.bam.bai
done 