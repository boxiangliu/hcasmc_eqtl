# WASP find_intersecting_snps:
python WASP_find_intersecting_snps.py ../data/eckersley/alignment/ sample_list.txt
python $scripts/WASP_find_intersecting_snps.py \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments/ \
	sample_list.txt


# WASP remap:
python WASP_remap.py ../data/eckersley/alignment sample_list_no_ERR440421.txt 


# WASP filter_remapped_reads: 
python WASP_filter_remapped_reads.py ../data/eckersley/alignment sample_list_no_ERR440421.txt 


# WASP rmdup: 
python WASP_rmdup.py ../data/eckersley/alignment sample_list_no_ERR440421.txt 


# add read group: 
python add_read_group.py ../data/eckersley/alignment sample_list_no_ERR440421.txt 
