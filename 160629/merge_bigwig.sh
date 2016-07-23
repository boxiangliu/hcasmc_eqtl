# merge bigwig to bedGraph: 
/software/ucsc_tools/3.0.9/bigWigMerge \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/bigwig/*/Aligned.out.sorted.rg.uniq.bw \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/bigwig/merged.bedGraph 

# convert bedGraph to bigwig: 
/software/ucsc_tools/3.0.9/bedGraphToBigWig \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/bigwig/merged.bedGraph  \
	/srv/persistent/bliu2/shared/genomes/hg19.chrom.sizes \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/bigwig/merged.bw
