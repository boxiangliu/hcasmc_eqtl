in_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/colocalization/HCASMC.clpp1e-2.txt.bak'
ecaviar_output_dir='../processed_data/eCAVIAR/eCAVIAR_output4/HCASMC'


# library
library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)

# read colocalized loci:
coloc=fread(in_file)
setnames(coloc,c('sid','score','rsid'))
coloc


# find gene name: 
genes=c()
for (i in 1:nrow(coloc)){
	sid=coloc[i,sid]
	score=coloc[i,score]
	stdout=system(sprintf('grep %s %s/*col',sid,ecaviar_output_dir),intern=TRUE)
	stdout_bak=stdout
	stdout=str_split_fixed(stdout_bak,'\t',n=2)
	possible_files=stdout[,1]
	coloc_scores=as.numeric(stdout[,2])
	coloc_file=possible_files[coloc_scores==score]
	tmp=str_split_fixed(basename(coloc_file),'.ecaviar_col:',n=2)
	gene=tmp[,1]
	genes=c(genes,gene)
}
coloc$gene=genes


# write output: 
write.table(coloc,'../processed_data/eCAVIAR/locuszoom/HCASMC.txt',quote=F,sep='\t',row.names=F,col.names=F)



