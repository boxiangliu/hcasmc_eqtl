library(XLConnect)


##############
# string to upper case function
###############

CapLeading <- function (string){
   fn <- function(x){
     v <- unlist(strsplit(x, split = " "))
     u <- sapply(v, function(x){
            x <- tolower(x)
            substring(x, 1, 1) <- toupper(substring(x, 1, 1))
            x})
     paste(u, collapse = " ")
     }
   sapply(string, fn)
   }

renameEthnicity = function(x){
	x = sub('african_american', 'AFR', x, fixed = TRUE)
	x = sub('caucasian', 'EUR', x, fixed = TRUE)
	x = sub('asian', 'EAS', x, fixed = TRUE)
	x = sub('hispanic', 'AMR', x, fixed = TRUE)
	x[which(is.na(x))] = 'unknown'
	return(x)
}


#---- main -----
# read command line args: 
args = commandArgs(T)
QFilename = args[1]
famFileName = args[2]
figure_name = args[3]

# read Q file (the ancestry proportion file):
# QFilename = '../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.biallelic.nomissing.noSAS.4.Q' 
QFile<-read.csv(QFilename,sep=" ",header=FALSE)
stopifnot(ncol(QFile) == 4)

# read fam file (output from plink): 
# famFileName = '../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.biallelic.nomissing.noSAS.fam'
famFile = read.csv(famFileName,sep=" ",header=FALSE)
stopifnot(nrow(famFile) == nrow(QFile))

# read sample sheet: 
hcasmc_panel = readWorksheet(loadWorkbook("../processed_data/R21R33_Cell_line_summary.xlsx"),sheet=1)
hcasmc_panel = as.data.table(hcasmc_panel)
setnames(hcasmc_panel, c('Catalog..', 'X..vials', 'gDNA.','RNA.'),c('Catalog', 'Vials', 'gDNA','RNA'))

# read 1000G ethnicity: 
TKG_panel = fread('../processed_data/integrated_call_samples_v3.20130502.ALL.panel', header = TRUE)
setnames(TKG_panel, c('sample','super_pop'), c('ID', 'Ethnicity'))


# row bind TGK panel and sample file: 
hcasmc_TGK_ethnicity = rbind(hcasmc_panel[,.(ID, Ethnicity)],TKG_panel[,.(ID, Ethnicity)])
stopifnot(nrow(hcasmc_TGK_ethnicity) == (2504 + 67 - 1))
hcasmc_TGK_ethnicity$Ethnicity2 = renameEthnicity(hcasmc_TGK_ethnicity$Ethnicity)

# column bind Q file and fam file: 
QFileFam = QFile
QFileFam = cbind(data.frame(ID = famFile[,1],QFile)) # this assumes that the fam file and Q file are in the same order

# merge QFileFam with ethnicity:  
QFileFamEthnicity<-merge(QFileFam,hcasmc_TGK_ethnicity, by = 'ID')


# change columns to population names: 
which(QFileFamEthnicity$Ethnicity2 == 'EUR')[1] # 1
QFileFamEthnicity[1,]
#        ID       V1    V2      V3       V4 Ethnicity Ethnicity2
# 1 1020301 0.027273 1e-05 0.01091 0.961808 caucasian        EUR

which(QFileFamEthnicity$Ethnicity2 == 'AFR')[1] # 15
QFileFamEthnicity[15,] 
#      ID       V1    V2    V3       V4        Ethnicity Ethnicity2
# 15 1497 0.872919 1e-05 1e-05 0.127061 african_american        AFR

which(QFileFamEthnicity$Ethnicity2 == 'EAS')[2] # 43
QFileFamEthnicity[45,] 
#      ID       V1       V2       V3       V4 Ethnicity Ethnicity2
# 45 2477 0.004393 0.921612 0.049447 0.024547     asian        EAS

which(QFileFamEthnicity$Ethnicity2 == 'AMR')[3] # 38
QFileFamEthnicity[38,] 
#      ID       V1       V2       V3       V4 Ethnicity Ethnicity2
# 38 2282 0.065446 0.009829 0.574028 0.350697  hispanic        AMR


population_columns = c('AFR','EAS','AMR','EUR')
setnames(QFileFamEthnicity, paste0('V', seq=1:4), population_columns)
KPlot<-t(QFileFamEthnicity[,population_columns])


## plot ancestry proportions (using bar plots): 
KNum = 4 # number of source populations. 
num_hcasmc = nrow(hcasmc_panel)
names = paste(QFileFamEthnicity$ID, QFileFamEthnicity$Ethnicity2, sep = '-')
pdf(figure_name, width = 10)
par(mai=c(1.5,0,.1,0.25),font=2,font.lab=4)
barplot(KPlot[,1:num_hcasmc], xlim = c(0,85), names.arg = names[1:num_hcasmc], col=rainbow(KNum),cex.names=.75, las = 2)
legend('right', col = rainbow(KNum), legend = population_columns, pch = 15)
dev.off()


