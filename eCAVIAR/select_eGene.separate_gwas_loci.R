# library: 
library(data.table)
library(dplyr)


# commandargs: 
args=commandArgs(T)
in_file=args[1]
out_dir=args[2]
chr=args[3]
# in_file='../processed_data/160816/subsampling/Liver/Liver_52.allpairs.sid_parsed.20.txt'
# out_dir='../processed_data/eCAVIAR/eCAVIAR_input/Liver/'
# chr=20

# read eqtl:
message(in_file)
eqtl=fread(in_file)
setnames(eqtl,c('pid','chr','pos','ref','alt','dist','pval_eqtl','beta_eqtl','se_eqtl'))


# if chr column of type character, coerce chr column of gwas to be character: 
if (typeof(eqtl$chr)=='character') gwas$chr=as.character(gwas$chr)


# read input file list:
gwas_files=list.files('../processed_data/eCAVIAR/cad_gwas_loci/',pattern=sprintf("^%s_",chr),full.names=T)
message(length(gwas_files),' gwas loci on chr',chr)

# read gwas:
for (gwas_file in gwas_files){
	gwas=fread(gwas_file)
	gwas=gwas[,.(chr, bp_hg19, markername, effect_allele, noneffect_allele, beta, se_dgc, p_dgc)]
	setnames(gwas,c('bp_hg19','beta','se_dgc','p_dgc'),c('pos','beta_gwas','se_gwas','pval_gwas'))


	# merge eqtl and gwas data: 
	merge=merge(eqtl,gwas,by=c('chr','pos'))
	merge[,ref2:=ifelse(nchar(ref)==1&nchar(alt)==1,ref,'I')]
	merge[,alt2:=ifelse(nchar(ref)==1&nchar(alt)==1,alt,'D')]
	merge=merge[pmin(ref2,alt2)==pmin(effect_allele,noneffect_allele)&pmax(ref2,alt2)==pmax(effect_allele,noneffect_allele)]


	# report: 
	n_total=length(unique(merge[,markername]))
	n_excluded=length(unique(merge[!(pmin(ref2,alt2)==pmin(effect_allele,noneffect_allele)&pmax(ref2,alt2)==pmax(effect_allele,noneffect_allele)),markername]))
	message(n_excluded,'/',n_total,' excluded because gwas and eqtl genotype differ!')
	# notes for hcasmc:
	# length(unique(merge[!(pmin(ref2,alt2)==pmin(effect_allele,noneffect_allele)&pmax(ref2,alt2)==pmax(effect_allele,noneffect_allele)),markername])) # 20
	# length(unique(merge[,markername])) # 13381
	# 20/13381 variants were exclued


	# select eGenes at GWAS loci (pval < 1e-3):
	selected_pid=unique(merge[pval_eqtl<1e-3,pid])
	message(length(selected_pid),' eGenes with p-value < 1e-3')
	selected=merge[pid%in%selected_pid,]


	# calculate zscore: 
	selected[,zscore_eqtl:=beta_eqtl/se_eqtl]
	selected[,zscore_gwas:=beta_gwas/se_gwas]


	# add snp id: 
	selected[,sid:=paste(chr,pos,ref,alt,'b37',sep='_')]


	# write output: 
	if (!dir.exists(out_dir)) dir.create(out_dir)
	for (id in selected_pid){
		tmp=selected[pid==id]
		tmp=tmp%>%arrange(chr,pos)%>%as.data.table()
		write.table(tmp[,.(sid,zscore_eqtl)],sprintf('%s/%s.eqtl.zscore',out_dir,id),quote=F,row.names=F,col.names=F,sep='\t')
		write.table(tmp[,.(sid,zscore_gwas)],sprintf('%s/%s.gwas.zscore',out_dir,id),quote=F,row.names=F,col.names=F,sep='\t')
		write.table(tmp[,sid],sprintf('%s/%s.id',out_dir,id),quote=F,row.names=F,col.names=F,sep='\t')
	}
}