library(data.table)
library(cowplot)
library(stringr)

in_dir='../processed_data/hcasmc_specific_open_chromatin/peak_specificity_filt/'
in_fn_ls=list.files(in_dir,pattern='bed')

pdf('../figures/hcasmc_specific_open_chromatin/distribution_of_tissue_sharing.pdf')
for (in_fn in in_fn_ls){
	sample=str_replace(in_fn,'.bed','')
	print(sample)
	peak=fread(sprintf("%s/%s",in_dir,in_fn))
	p=ggplot(peak,aes(num_tissue))+geom_histogram(binwidth=1)+ggtitle(sample)
	print(p)
}
dev.off()