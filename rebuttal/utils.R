library(data.table)
library(XLConnect)

covariate_fn = '../data/sample_info/sample_info.xlsx'
read_known_covariate = function(covariate_fn){
	covariate = readWorksheet(loadWorkbook(covariate_fn),sheet=5)
	setDT(covariate)
	return(covariate)
}
