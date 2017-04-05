library(data.table)
library(stringr)

# Variables: 
in_dir='../processed_data/hcasmc_specific_open_chromatin/peak_specificity/'
out_dir='../processed_data/hcasmc_specific_open_chromatin/peak_specificity_filt/'
if (!dir.exists(out_dir)) dir.create(out_dir)

# Get file names: 
fn_list=list.files(in_dir,'.bed')


# Loop through each file: 
for (fn in fn_list){
	print(fn)
	peak=fread(sprintf('%s/%s',in_dir,fn),header=F,select=c(1:4,8,9),col.names=c('chr','start','end','signalValue','id','psi'))
	
	peak[,num_tissue:=as.integer(str_split_fixed(id,'_',4)[,4])]
	peak[,id:=NULL]

	peak[,min_psi:=min(psi),by=c('chr','start','end')]
	peak2=unique(peak[psi==min_psi,])
	
	peak2[,min_psi:=NULL]

	out_fn=sprintf('%s/%s',out_dir,fn)
	fwrite(peak2,out_fn,sep='\t')
}