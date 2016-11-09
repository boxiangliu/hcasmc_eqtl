# extract GWAS hits from the supplmentary excel in the paper Nikpay 2015 Nature Genetics
library(XLConnect)
library(jsonlite)


# command args: 
args=commandArgs(T)
in_file=args[1]
out_file=args[2]
# in_file="/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/nikpay_2015_ng.xlsx"
# out_file='../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt'

# functions: 
rsid2pos=function(rsids){
	response=data.frame()
	for (rsid in rsids){
		tmp=fromJSON(system(sprintf("curl 'http://grch37.rest.ensembl.org/variation/human/%s?' -H 'Content-type:application/json'",rsid),intern=T))$mappings[1,]%>%select(chr=seq_region_name,pos=start)
		response=rbind(response,tmp)
	}
	return(response)
}


# read known gwas loci:
gwas_known_hits=readWorksheet(loadWorkbook(in_file),sheet=3,startRow=3,endRow=58,startCol=2,endCol=3,check.names=F,header=T)%>%select(chr=Chr,markername=`Published SNP`)


# get position for gwas known hits:
gwas_known_hits$pos=rsid2pos(gwas_known_hits$markername)$pos


# read FDR hits: 
gwas_fdr_hits=readWorksheet(loadWorkbook(in_file),sheet=7,startRow=2,endRow=204,startCol=1,endCol=3,check.names=F,header=T)%>%select(markername,chr,pos=bp_hg19)


# read genome wide significant hits (CAD): 
gwas_genowide_hits_cad=readWorksheet(loadWorkbook(in_file),sheet=6,startRow=3,endRow=59,startCol=2,endCol=4,check.names=F,header=T)%>%select(markername=SNP,chr=CHR,pos=`Base-pair Position`)


# read genome wide significant hits (MI):
gwas_genowide_hits_mi=readWorksheet(loadWorkbook(in_file),sheet=6,startRow=3,endRow=59,startCol=2,endCol=11,check.names=T,header=T)%>%select(markername= SNP.1,chr=CHR,pos=`Base.pair.Position.1`)


# read genome wide lead hits (both additive and recessive): 
gwas_lead_variants=readWorksheet(loadWorkbook(in_file),sheet=5,startRow=4,endRow=14,startCol=1,endCol=3,check.names=F,header=T)%>%select(markername=`Lead variant`,chr=Chr.)
gwas_lead_variants$pos=rsid2pos(gwas_lead_variants$markername)$pos


# merge two lists:
gwas_known_hits=gwas_known_hits%>%select(markername,chr,pos)
gwas_fdr_hits=gwas_fdr_hits%>%select(markername,chr,pos)
gwas_genowide_hits_cad=gwas_genowide_hits_cad%>%select(markername,chr,pos)
gwas_genowide_hits_mi=gwas_genowide_hits_mi%>%select(markername,chr,pos)
gwas_lead_variants=gwas_lead_variants%>%select(markername,chr,pos)
merged=rbind(gwas_known_hits,gwas_fdr_hits,gwas_genowide_hits_cad,gwas_genowide_hits_mi,gwas_lead_variants)


# take unique hits: 
uniq=unique(merged)

# output the merged unique hits:
write.table(uniq,out_file,sep='\t',row.names=F,quote=F)

	