# command args: 
args=commandArgs(T)
in_file=args[1]
figure_prefix=args[2]
out_file=args[3]

# in_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output/metasoft_output.5.mcmc.txt'
# figure_prefix='../figures/hcasmc_specific_eqtl/'
# out_file='../processed_data/hcasmc_specific_eqtl/hcasmc_specific_eqtl.chr5.txt'


# constants: 
N_TISSUE=45


# read metasoft result: 
header=scan(in_file,what='character',nlines=1,sep='\t')
metasoft=fread(in_file,skip=1)


# remove the extra column (issue due to white space)
n_col=ncol(metasoft) 
metasoft=metasoft[,1:(n_col-1),with=F]


# read tissue name (study name):
study_name=unlist(fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.alphabetical.txt',header=F))


# set metasoft column names:
col_names=c(paste('pvalue',study_name,sep='_'),paste('mvalue',study_name,sep='_'))
col_names=c(header[1:16],col_names)
stopifnot(ncol(metasoft)==length(col_names))
setnames(metasoft,col_names)


# subset to mvalues: 
mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
stopifnot(ncol(mvalue)==N_TISSUE)
rownames(mvalue)=metasoft$RSID
colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')



# subset to pvalues: 
pvalue=metasoft%>%dplyr::select(contains('pvalue'))%>%dplyr::select(-c(1:5))%>%as.data.frame()
rownames(pvalue)=metasoft$RSID
colnames(pvalue)=str_replace(colnames(mvalue),'pvalue_','')


# select rows where HCASMC has maximum m-value: 
hcasmc_highest=mvalue$HCASMC==apply(mvalue,1,max,na.rm=T)


# select rows where HCASMC has m-value >= 0.9:
hcasmc_eqtl=mvalue$HCASMC>=0.9


# select rows where GTEx tissues all have m-values < 0.9:
gtex_max=apply(mvalue%>%select(-HCASMC),1,max,na.rm=T)
gtex_non_eqtl=gtex_max<0.9


# select hcasmc-specific eQTLs:
hcasmc_specific_eqtl_idx=which(hcasmc_eqtl&hcasmc_highest&gtex_non_eqtl)
message(paste(length(hcasmc_specific_eqtl_idx),'hcasmc specific eqtls.'))

# make PM plot:
for (idx in hcasmc_specific_eqtl_idx){
	to_plot=data.frame(mvalue=unlist(mvalue[idx,]),pvalue=-log10(unlist(pvalue[idx,])))
	to_plot$tissue=rownames(to_plot)
	p=ggplot(to_plot,aes(mvalue,pvalue,label=tissue))+geom_point()+geom_text(angle=45)+theme_bw()+geom_vline(xintercept=0.9,color='red',linetype='dashed')
	save_plot(paste0(figure_prefix,rownames(mvalue[idx,]),'.pdf'),p,base_width=6,base_height=6)
}


# write hcasmc-specific eQTL to text file: 
output=as.data.frame(str_split_fixed(names(hcasmc_specific_eqtl_idx),"_",n=2))
setnames(output,c('pid','sid'))
write.table(output,file=out_file,sep='\t',row.names=F,quote=F,col.names=F)