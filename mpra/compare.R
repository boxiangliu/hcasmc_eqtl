library(data.table)
library(dplyr)
library(cowplot)
naive_file='../processed_data/mpra/naive_overlap_eQTL/gwas_eQTL_naive_overlap.txt'
ecaviar_file='../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.allFeatures.bed'
naive=fread(naive_file)
ecaviar=fread(ecaviar_file)

cuts=c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
num_var=c()
for (cut in cuts){
	num_var=c(num_var,naive[pval<cut,rsid]%>%unique()%>%length())
}
ggplot(data.frame(cbind(cuts,num_var)),aes(-log10(cuts),num_var))+geom_point()+geom_text(label=num_var,nudge_y=50)+xlab('P-value cutoff')+ylab('Number of variants')
