# constants: 
N_TISSUE=45


# read metasoft result: 
header=scan('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output/metasoft_output.5.mcmc.txt',what='character',nlines=1,sep='\t')
metasoft=fread('/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_output/metasoft_output.5.mcmc.txt',skip=1)


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


metasoft[which(metasoft$RSID=='ENSG00000197208.5_5_131667353_A_G_b37'),]
# subset to mvalues: 
mvalue=metasoft%>%dplyr::select(contains('mvalue'))%>%as.data.frame()
stopifnot(ncol(mvalue)==N_TISSUE)


# set row and column names for mvalues: 
rownames(mvalue)=metasoft$RSID
colnames(mvalue)=str_replace(colnames(mvalue),'mvalue_','')


# subset to rows where HCASMC mvalue > 0.9: 
hcasmc_eqtl_idx=which(mvalue$HCASMC>0.9)
max_mvalue=apply(hcasmc%>%select(-HCASMC),1,max,na.rm=T)
hcasmc[hcasmc$HCASMC>max_mvalue,]
which(hcasmc$HCASMC>max_mvalue)


# subset to pvalues: 
pvalue=metasoft%>%dplyr::select(contains('pvalue'))%>%dplyr::select(-c(1:5))%>%as.data.frame()


# set row and column names for mvalues: 
rownames(pvalue)=metasoft$RSID
colnames(pvalue)=str_replace(colnames(mvalue),'pvalue_','')


# PM plot:

idx=hcasmc_eqtl_idx[1]
to_plot=data.frame(mvalue=unlist(mvalue[idx,]),pvalue=-log10(unlist(pvalue[idx,])))
to_plot$tissue=rownames(to_plot)
ggplot(to_plot,aes(mvalue,pvalue,label=tissue))+geom_point()+geom_text(angle=45)

to_plot