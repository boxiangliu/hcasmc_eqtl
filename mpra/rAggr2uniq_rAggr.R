library(data.table)
library(dplyr)
library(stringr)
args=commandArgs(T)
raggr_file=args[1]
raggr_output_file=args[2]
rsq=args[3]

# raggr_file='../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.6_1kgp3.csv'
# raggr_output_file='../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.8_1kgp3.uniq.txt'
# rsq=0.8

raggr=fread(raggr_file)%>%select(chr=`SNP2 Chr`,name=`SNP2 Name`,r2=`R-squared`)
raggr_rsq=raggr%>%filter(r2>=rsq)%>%select(chr,name)
raggr_uniq=unique(raggr_rsq)
split_name=str_split_fixed(raggr_uniq$name,':',n=4)
raggr_uniq=cbind(raggr_uniq,split_name)
setnames(raggr_uniq,3:6,c('markername_bak','pos','ref','alt'))
raggr_output=raggr_uniq%>%mutate(markername=ifelse(str_detect(markername_bak,'rs'),as.character(markername_bak),paste(chr,pos,ref,alt,sep=':')))%>%select(markername,chr,pos)
write.table(raggr_output,raggr_output_file,row.names=F,sep='\t',quote=F)