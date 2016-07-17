#!/usr/bin/Rscript
# 6/1/16
# Joe Davis
# Modified by B Liu
# R script to choose which individuals to subsample for in each tissue at each sample size

# Define subsample sizes
args=commandArgs(T,T)
subsample.size=args$size
subsample.size=10

# Read in the sample annotations and subset to 
# samples with read count files: 
annotation=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt',sep='\t',header=T,stringsAsFactors=F)
samples_with_read_count=t(read.table('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_All_Tissues_read_counts_FOR_QC_ONLY.gct',skip=2,nrows=1)[,-c(1,2)])
stopifnot(samples_with_read_count%in%annotation$SAMPID)
stopifnot(unique(annotation$SAMPID)==annotation$SAMPID)
annotation=annotation[SAMPID%in%samples_with_read_count,]


# reformat the SMTSD column (tissue type) of annotation data.frame: 
annotation[,SMTSD2:=SMTSD%>%str_replace_all("- |\\(|\\)","")%>%str_replace_all(" |-",'_')%>%str_replace_all(" ","")]
annotation[,NSAMP:=.N,by='SMTSD2']
annotation=annotation[NSAMP>=subsample.size,]
tissues=unique(annotation$SMTSD2)


# For each tissue, read in the list of individuals and subsample that list down until the smallest subsample size is reached
# Output each list to the correct tissue directory for use in the FastQTL pipeline
set.seed(1)
for(i in 1:length(tissues)){
	message(tissues[i])
	tissue.samples = annotation[SMTSD2==tissues[i],SAMPID]
	message(length(tissue.samples))
	subsample=sort(sample(tissue.samples, subsample.size))
	output_dir=paste0('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/', 'subsampling/', tissues[i])
	if (!dir.exists(output_dir)){dir.create(output_dir)}
	write.table(subsample, paste(output_dir, '/', tissues[i], '_subsample.', subsample.size, '.txt', sep = ''), col.names = F, row.names = F, quote = F)
}

