library(XLConnect)
library(data.table)

# Correlate PEER covariates with age and sex: 
covariate_fn = '../data/sample_info/sample_info.xlsx'
peer_fn = '../processed_data/eqtl/peer/factors.tsv'

read_known_covariate = function(covariate_fn){
	covariate=readWorksheet(loadWorkbook(covariate_fn),sheet=5)
	setDT(covariate)
	return(covariate)
}

covariate = read_known_covariate(covariate_fn)

read_peer_covariate = function(peer_fn){
	peer = read.table(peer_fn,header=T,row.names=1,check.names=F)
	peer = data.frame(t(peer))
	peer$sample = rownames(peer)
	setDT(peer)
	return(peer)
}

peer = read_peer_covariate(peer_fn)
merged = merge(peer,covariate[,list(RNA = as.character(RNA),Sex,Age)],by.x='sample',by.y='RNA')

merged[,Sex:=ifelse(Sex=='M',0,1)]
correlation = cor(merged[,2:ncol(merged)])
max(abs(correlation['Age',-17])) # 0.3565919
max(abs(correlation['Sex',-16])) # 0.3397523