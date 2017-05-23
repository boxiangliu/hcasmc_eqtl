library(XLConnect)
library(flashpcaR)
library(gplots)
source('utils.R')



renameEthnicity = function(x){
	x = sub('african_american', 'AFR', x, fixed = TRUE)
	x = sub('caucasian', 'EUR', x, fixed = TRUE)
	x = sub('asian', 'EAS', x, fixed = TRUE)
	x = sub('hispanic', 'AMR', x, fixed = TRUE)
	x[which(is.na(x))] = 'unknown'
	return(x)
}


#------ main ------
args = commandArgs(TRUE)
inp = args[1]
fig = args[2]


# read genotypes: 
inp = "../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.vcf.GT.FORMAT"
genotypes = fread(inp, header = T)
genotypes[,'150328' := NULL,with=F]
genotypes = removeMissingGenotypes(genotypes)
genotypes = removeMultiAllelicLoci(genotypes) # remove multi allelic loci
dosages = genotypeToDosage(genotypes) # convert genotypes to dosages



# read sample sheet: 
hcasmc_panel = readWorksheet(loadWorkbook("../processed_data/R21R33_Cell_line_summary.xlsx"),sheet=1)
hcasmc_panel = as.data.table(hcasmc_panel)
setnames(hcasmc_panel, c('Catalog..', 'X..vials', 'gDNA.','RNA.'),c('Catalog', 'Vials', 'gDNA','RNA'))

# read 1000G ethnicity: 
TKG_panel = fread('../processed_data/integrated_call_samples_v3.20130502.ALL.panel', header = TRUE)
setnames(TKG_panel, c('sample','super_pop'), c('ID', 'Ethnicity'))

# merge hcasmc sample sheet and 1000G panel:  
colData = rbind(hcasmc_panel[,.(ID, Ethnicity)],TKG_panel[,.(ID, Ethnicity)])
stopifnot(nrow(colData) == (2504 + 67 - 1))

# order colData to match dosages column names: 
idx = match(colnames(dosages), colData$ID)
colData = colData[idx,]
stopifnot(colData$ID == colnames(dosages))
stopifnot(unique(colData$ID) == colData$ID)

# remove South Asian: 
idx = which(colData$Ethnicity == 'SAS')
colData = colData[-idx,]
dosages = dosages[,-idx]
stopifnot(colData$ID == colnames(dosages))

# remove Asian:
idx = which(colData$Ethnicity == 'EAS')
colData = colData[-idx,]
dosages = dosages[,-idx]
stopifnot(colData$ID == colnames(dosages))

# plot PCs: 
r = flashpca(dosages, do_loadings=TRUE, verbose=TRUE, stand="binom", ndim=10,nextra=100)
pcs = r$loadings
colnames(pcs) = paste0('PC', seq(1:10))
pcs = data.table(pcs, Ethnicity = colData$Ethnicity)
pcs[, size := ifelse(Ethnicity %in% c('AFR','EUR','EAS','SAS','AMR'), 10, 20), by = Ethnicity]
pcs[, alpha := ifelse(Ethnicity %in% c('AFR','EUR','EAS','SAS','AMR'), 0.0, 1), by = Ethnicity]
pcs[, color := ifelse(Ethnicity %in% c('AFR','EUR','EAS','SAS','AMR'), 0.1, 1), by = Ethnicity]
pcs[, Ethnicity2 := Ethnicity ]
Ethnicity2 = renameEthnicity(colData$Ethnicity)
pcs$Ethnicity2 = Ethnicity2

ggplot(data = pcs, aes(PC1, PC2, color = Ethnicity2, alpha = alpha)) + geom_point()
ggplot(data = pcs, aes(PC2, PC3, color = Ethnicity2, alpha = alpha)) + geom_point()
ggplot(data = pcs, aes(PC1, PC3, color = Ethnicity2, alpha = alpha)) + geom_point()

# load saved image: 
save.image('../r_environments/plot_ancestry_PCs.RData')
load('../r_environments/plot_ancestry_PCs.RData')