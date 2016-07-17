#!/usr/bin/env Rscript
# bosh liu
# durga
# check that the ref and alt alleles are on the same strand for two input files

# command args: 
args=commandArgs(T)
input1_file=args[1]
input2_file=args[2]
input1_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603_phasing/chr22.1kg.phase3.v5a.txt'
input2_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603_phasing/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.GRCh37.txt'

# read inputs: 
input1=fread(input1_file,header=T)
input2=fread(input2_file,header=T)


# set column names: 
setnames(input1,c('CHROM','POS','REF','ALT'))
setnames(input2,c('CHROM','POS','REF','ALT'))


# subset to chr22:
input2=input2[CHROM=='22',]


# make ID column: 
input1[,ID:=paste(CHROM,POS,sep="_")]
input2[,ID:=paste(CHROM,POS,sep="_")]


# merge input1 and input2:
merged=merge(input1[,.(ID,REF,ALT)],input2[,.(ID,REF,ALT)],by='ID')
merged[,same_ref:=(REF.x==REF.y)]
merged[,same_alt:=(ALT.x==ALT.y)]
merged[same_ref==F,.(REF.x,REF.y)]
merged[same_alt==F,.(ID,REF.x,REF.y,ALT.x,ALT.y)]

