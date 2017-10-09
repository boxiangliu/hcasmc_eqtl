# Variables:
bed_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_released_all_cell_type/'
out_dir=''


# Compile file list: 
bed_fn_list=list.files(bed_dir,pattern='bed',full.names=TRUE)

# Calculate jaccard similarity using Bedtools: 
i=1
j=1
fa=bed_fn_list[i]
fb=bed_fn_list[j]
cmd=sprintf('bedtools jaccard -a %s -b %s',fa,fb)

# Remove tissue with high jaccard similarity:


