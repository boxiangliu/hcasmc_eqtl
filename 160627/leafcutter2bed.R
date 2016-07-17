# commandargs:
args=commandArgs(T,T)
leafcutter_file=args$leafcutter
bed_file=args$bed
leafcutter_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.tsv'
bed_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.bed'


# read splice levels:
leafcutter=fread(leafcutter_file,header=T)


# parse chromosome and location:
loci=str_split_fixed(leafcutter$Name,":",4)[,1:3]
loci=as.data.table(loci)
setnames(loci,c("#chr","start","end"))


# generate bed file:
bed=data.table(loci,leafcutter)
setnames(bed,'Name','ID')


# sort by chr, start and end:
bed$chr=str_replace(bed$`#chr`,'chr','')
bed[chr=='X',]$chr='23'
bed[chr=='Y',]$chr='24'
bed=arrange(bed,as.numeric(chr),as.numeric(start),as.numeric(end))
bed[,chr:=NULL]


# output bed file:
write.table(bed,bed_file,quote=F,sep='\t',row.names=F,col.names=T)