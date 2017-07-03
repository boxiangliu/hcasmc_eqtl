rm /srv/persistent/bliu2/HCASMC_eQTL/scripts//atacseq/wasp/sample_file.txt 
for s in 1346 1508 1522 200212 2108 2305 2356 2510 2989; do
	echo $s
	fastq1=`ls /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/$s/out/align/rep1/*R1*.fastq.gz`
	fastq2=`ls /srv/persistent/bliu2/HCASMC_eQTL/data/atacseq/fbs/$s/out/align/rep1/*R2*.fastq.gz`
	echo $s $s $fastq1 $fastq2 >> /srv/persistent/bliu2/HCASMC_eQTL/scripts//atacseq/wasp/sample_file.txt 
done 