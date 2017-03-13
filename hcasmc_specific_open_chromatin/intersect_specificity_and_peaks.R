# Variables: 
peak_dir='../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/'
spec_dir='../processed_data/hcasmc_specific_open_chromatin/specificity_index/'
out_dir='../processed_data/hcasmc_specific_open_chromatin/peak_specificity/'

if (!dir.exists(out_dir)){dir.create(out_dir)}

# Get all peak files: 
peak_fn=list.files(peak_dir,pattern='bed')
spec_fn=list.files(spec_dir,pattern='bed')
stopifnot(all(peak_fn==spec_fn))


# Loop through all samples: 
for (f in peak_fn){
	# Run bedtools: 
	print(f)
	command=sprintf("bedtools intersect -a %s/%s -b %s/%s -wa -wb > %s/%s",peak_dir,f,spec_dir,f,out_dir,f)
	system(command)
}