## library:
library(data.table)
library(qvalue)
library(dplyr)
library(VennDiagram)

## variables: 
fastqtl_file='../processed_data/compare_rasqual_and_fastqtl/fastqtl.txt'
rasqual_file='../processed_data/compare_rasqual_and_fastqtl/rasqual.txt'
out_fig='../figures/compare_rasqual_and_fastqtl/venn.pdf'
if (!dir.exists(dirname(out_fig))) {dir.create(dirname(out_fig),recursive=T)}


## main: 
# Read files:
fastqtl=fread(fastqtl_file)
rasqual=fread(rasqual_file)
setnames(rasqual,c('fid','chr_pos_ref_alt'),c('gene_name','sid'))
merged=merge(fastqtl,rasqual,by=c('gene_name','sid'))
setnames(merged,c('pval.x','qval.x','pval.y','qval.y'),c('pval_fastqtl','qval_fastqtl','pval_rasqual','qval_rasqual'))


# Recalculate Q-values to account for dropped entries during merging: 
merged[,qval_fastqtl:=NULL]
merged[,qval_rasqual:=NULL]
merged$qval_fastqtl=qvalue(merged$pval_fastqtl)$qvalues
merged$qval_rasqual=qvalue(merged$pval_rasqual)$qvalues


# Calculate the number of discoveries using fastQTL and RASQUAL: 
fq_dis=sum(merged$qval_fastqtl<0.05)
rq_dis=sum(merged$qval_rasqual<0.05)
joint_dis=sum(merged$qval_fastqtl<0.05&merged$qval_rasqual<0.05)
rq_dis/fq_dis
joint_dis/fq_dis


# Make Venn plot to compare RASQUAL and fastQTL discoveries: 
pdf(out_fig)
draw.pairwise.venn(fq_dis,rq_dis,joint_dis,category = c("Total reads only", "Total reads and ASE"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()


# Calculate the proportion of variants with same directory of effect: 
joint=merged[qval_fastqtl<0.05&qval_rasqual<0.05,]
joint[,concord:=!xor(pi>0.5,beta>0)]
nrow(joint[concord==1,])/nrow(joint)

