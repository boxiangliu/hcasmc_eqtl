# bosh liu
# 16/08/06
# make metasoft input file: 

# library: 
library(qvalue)


# read in tissue name: 
tissue_filename="/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/Metasoft_tissue_order.txt"
tissues=fread(tissue_filename,header=F)%>%unlist()%>%unname()
filename_format='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/%s_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz'
# filename_format='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/%s_Analysis.v6p.FOR_QC_ONLY.allpairs.head2000.txt'


master=data.table()
for (tissue in tissues){
	# report progress: 
	message(tissue)

	# read eqtl data:
	eqtl_file=sprintf(filename_format,tissue)
	if (str_detect(eqtl_file,'\\.gz')) {
		eqtl=fread(sprintf('zcat %s',eqtl_file))
	} else {
		eqtl=fread(eqtl_file)
	}

	# set names: 
	setnames(eqtl,c('pheno','geno','dist','pval','beta','se'))


	# remove dist column: 
	eqtl[,dist:=NULL]


	# subset to one gene for testing: 
	# eqtl=eqtl[pheno=='ENSG00000227232.4',]


	# add ID column: 
	eqtl[,ID:=paste(pheno,geno,sep='_')]


	# remove pheno and geno column: 
	eqtl[,pheno:=NULL]
	eqtl[,geno:=NULL]


	# calcuale qvalue: 
	# eqtl$qval=qvalue(eqtl$pval)$qvalues


	# remove pvalue: 
	# eqtl[,pval:=NULL]


	# add tissue name to beta, se and qval column: 
	# setnames(eqtl,c('beta','se','qval'),paste(c('beta','se','qval'),tissue,sep="_"))
	setnames(eqtl,c('beta','se','pval'),paste(c('beta','se','pval'),tissue,sep="_"))


	if (ncol(master)==0) {
		master=eqtl
	} else {
		master=merge(master,eqtl,by='ID',all=TRUE)
	}
}



# select pairs with at least tissue significant at qvalue 5%: 
# qval=master%>%dplyr::select(contains('qval'))
# min_qval=apply(qval,1,min,na.rm=TRUE)
# selected=which(min_qval<0.05)
# master=master[selected,]

pval=master%>%dplyr::select(contains('pval'))
min_pval=apply(pval,1,min,na.rm=TRUE)
selected=which(min_pval<1e-3)
master=master[selected,]

# remove qvalues columns: 
# output=master%>%dplyr::select(-contains('qval'))
output=master%>%dplyr::select(-contains('pval'))


# write output:
write.table(output,'/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/metasoft_input.txt',quote=F,col.names=F,row.names=F)
