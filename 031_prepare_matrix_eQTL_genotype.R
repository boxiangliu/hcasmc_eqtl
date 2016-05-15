#!/usr/bin/env Rscript
# bosh liu
# 2016/05/08
# durga
# prepare data in matrix eQTL format

# library:
library(XLConnect)
library(data.table)

# command args: 
args=commandArgs(T)
chr=args[1]


# path:
snp_file=sprintf('../processed_data/028_imputed_genotype/chr%s.genotype',chr)
message(snp_file)
output_dir='../processed_data/031_prepare_matrix_eQTL_genotype/'



# read sample_sheet:
sample_sheet=readWorksheet(loadWorkbook("../processed_data/rna_wgs_match.reduced_050616.xlsx"),sheet=1)
sample_sheet=as.data.table(sample_sheet)



# remove 9052004_dase, because 
# 9052005_dase and 9052004 were merged in expression data
sample_sheet=sample_sheet[RNA.New.Name!='9052004_dase']



# subset to useful columns:
sample_sheet=sample_sheet[!is.na(RNA.New.Name),.(dna,RNA.New.Name,DNA.New.Name)]



# read snp data
snp=fread(snp_file,header=T)



# subset snp columns to working set:
snp=snp[,colnames(snp)%in%sample_sheet[,dna]|colnames(snp)=='SNP',with=F]



# rearrange snp columns to match working set:
snp=snp[,c('SNP',sample_sheet[,dna]),with=F]


# rename snp columns to DNA.new.name: 
setnames(snp,as.character(sample_sheet[,dna]),as.character(sample_sheet[,DNA.New.Name]))


# sort snp columns alphanumerically: 
snp=snp[,c('SNP',sort(colnames(snp)[-1])),with=F]


# change column name from "SNP" to "id":
setnames(snp,'SNP','id')


# write table: 
if (!dir.exists(output_dir)){dir.create(output_dir)}
write.table(snp,file=sprintf('%s/chr%s.genotype.txt',output_dir,chr),row.names=F,quote=F,sep='\t')




