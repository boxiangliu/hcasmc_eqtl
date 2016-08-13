#########
# SETUP
#########
library(dplyr)
library(data.table)
library(DESeq2)
library("BiocParallel")
source('/srv/persistent/bliu2/HCASMC_eQTL/scripts/utils.R')
register(MulticoreParam(10))

#######
# MAIN
#######

####################
# read count matrix
####################
count_wide = fread('../processed_data/160801/rnaseq_dase.combined.count')

# remove the singleton 9090701: 
count_wide=count_wide%>%dplyr::select(-`9090701_sf`)

# coerce to data.frame: 
count_matrix = count_wide %>% dplyr::select(-Name,-Description) %>% data.frame(check.names=F)
rownames(count_matrix) = count_wide %>% select(Name) %>% unlist()


########################
# create colData matrix
########################
# individual:
individual = str_split_fixed(colnames(count_matrix),'_',2)[,1]

# condition: 
condition = str_split_fixed(colnames(count_matrix),'_',2)[,2] %>% str_replace(.,"_[123]",'')

# batch:
batch2 = c('9090701_sf','8072501_sf','8100901_sf','1051601_sf','9071501_sf','9052004_fbs','9052004_sf')
batch=ifelse(colnames(count_matrix) %in% batch2, 'batch2','batch1')
batch=ifelse(str_detect(colnames(count_matrix),'2305'),'batch3',batch)

# colData: 
colData = data.frame(
  row.names = colnames(count_matrix), 
  condition = condition,
  individual = individual,
  batch = batch
)


########################
# create DESeq2 dataset
########################
dds = DESeqDataSetFromMatrix(countData = count_matrix,
                             colData = colData,
                             design = ~ condition + batch)

# change base level to sf. 
dds$condition = relevel(dds$condition, 'sf') 

# change base level to batch1. 
# dds$batch = relevel(dds$batch, 'batch1')


##############
# run DESeq2 
##############
dds = DESeq(dds, parallel = TRUE)
res = results(dds, contrast = c('condition','fbs','sf'))
resOrdered = res[order(res$padj),]
sum(resOrdered$padj<1e-2,na.rm=T) # 2472
up=sum(with(resOrdered,padj<1e-2&log2FoldChange>0),na.rm=T) # 1380
down=sum(with(resOrdered,padj<1e-2&log2FoldChange<0),na.rm=T) # 1092
insig=sum(resOrdered$padj>1e-2,na.rm=T) # 20793
untest=sum(is.na(resOrdered$padj)) # 32973

# make MA plot: 
pdf('../figures/160801/HCASMC_SF_vs_FBS.pdf')
plotMA(res,alpha=1e-2,main='HCASMC log(FBS/SF)',ylim=c(-5,5))
text(1,4,sprintf('Up: %s\nDown: %s\nInsignifant: %s\nUntested: %s',up,down,insig,untest),col='black')
dev.off()


# inspect the top 10 genes: 
count_wide[count_wide$Name%in%rownames(resOrdered[1:10,]),.(Name,Description)]
#  1:  ENSG00000117399.9       CDC20
# CDC20 appears to act as a regulatory protein interacting with several other proteins at multiple points in the cell cycle. It is required for two microtubule-dependent processes, nuclear movement prior to anaphase and chromosome separation. [provided by RefSeq, Jul 2008]
#  2:  ENSG00000183856.6      IQGAP3
# IQGAP3 (IQ Motif Containing GTPase Activating Protein 3) is a Protein Coding gene. Among its related pathways are Signaling by GPCR and Signaling by Rho GTPases. GO annotations related to this gene include calmodulin binding and Ras GTPase binding. An important paralog of this gene is IQGAP1.
#  3:  ENSG00000117650.8        NEK2
# This gene encodes a serine/threonine-protein kinase that is involved in mitotic regulation. This protein is localized to the centrosome, and undetectable during G1 phase, but accumulates progressively throughout the S phase, reaching maximal levels in late G2 phase. Alternatively spliced transcript variants encoding different isoforms with distinct C-termini have been noted for this gene. [provided by RefSeq, Feb 2011]
#  4:  ENSG00000171848.9        RRM2
# This gene encodes one of two non-identical subunits for ribonucleotide reductase. This reductase catalyzes the formation of deoxyribonucleotides from ribonucleotides. Synthesis of the encoded protein (M2) is regulated in a cell-cycle dependent fashion. Transcription from this gene can initiate from alternative promoters, which results in two isoforms that differ in the lengths of their N-termini. Related pseudogenes have been identified on chromosomes 1 and X. [provided by RefSeq, Sep 2009]
#  5:  ENSG00000152253.4       SPC25
# This gene encodes a protein that may be involved in kinetochore-microtubule interaction and spindle checkpoint activity. [provided by RefSeq, Jul 2008]
#  6:  ENSG00000112984.7      KIF20A
# KIF20A (Kinesin Family Member 20A) is a Protein Coding gene. Diseases associated with KIF20A include charcot-marie-tooth disease, type 4c and charcot-marie-tooth disease. Among its related pathways are Platelet activation, signaling and aggregation and Immune System. GO annotations related to this gene include protein kinase binding and ATPase activity. An important paralog of this gene is KIF20B.
#  7:  ENSG00000126787.8      DLGAP5
# Potential cell cycle regulator that may play a role in carcinogenesis of cancer cells. Mitotic phosphoprotein regulated by the ubiquitin-proteasome pathway. Key regulator of adherens junction integrity and differentiation that may be involved in CDH1-mediated adhesion and signaling in epithelial cells.
#  8: ENSG00000089685.10       BIRC5
# This gene is a member of the inhibitor of apoptosis (IAP) gene family, which encode negative regulatory proteins that prevent apoptotic cell death. IAP family members usually contain multiple baculovirus IAP repeat (BIR) domains, but this gene encodes proteins with only a single BIR domain. The encoded proteins also lack a C-terminus RING finger domain. Gene expression is high during fetal development and in most tumors, yet low in adult tissues. Alternatively spliced transcript variants encoding distinct isoforms have been found for this gene. [provided by RefSeq, Jun 2011]
#  9: ENSG00000101057.11       MYBL2
# The protein encoded by this gene, a member of the MYB family of transcription factor genes, is a nuclear protein involved in cell cycle progression. The encoded protein is phosphorylated by cyclin A/cyclin-dependent kinase 2 during the S-phase of the cell cycle and possesses both activator and repressor activities. It has been shown to activate the cell division cycle 2, cyclin D1, and insulin-like growth factor-binding protein 5 genes. Two transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jul 2013]
# 10: ENSG00000175063.12       UBE2C
# The modification of proteins with ubiquitin is an important cellular mechanism for targeting abnormal or short-lived proteins for degradation. Ubiquitination involves at least three classes of enzymes: ubiquitin-activating enzymes, ubiquitin-conjugating enzymes, and ubiquitin-protein ligases. This gene encodes a member of the E2 ubiquitin-conjugating enzyme family. The encoded protein is required for the destruction of mitotic cyclins and for cell cycle progression, and may be involved in cancer progression. Multiple transcript variants encoding different isoforms have been found for this gene. Pseudogenes of this gene have been defined on chromosomes 4, 14, 15, 18, and 19. [provided by RefSeq, Aug 2013]


# rank genes in descending order based on pvalues: 
ranks = rank(-resOrdered$pvalue,ties.method='random',na.last=F)
temp=data.frame(gene_id=rownames(resOrdered),rank=ranks)
output=merge(temp,count_wide[,.(Name,Description)],by.x='gene_id',by.y='Name')
output=output%>%arrange(-rank)
output=output%>%dplyr::select(gene_name=Description,rank)



# output ranked list: 
write.table(output,file='../processed_data/160801/HCASMC_SF_vs_FBS.rnk',quote=F,row.names=F,sep='\t')


# save:
save(list = c('dds', 'res', 'resOrdered'), file = '../processed_data/160801/DESeq2_HCASMC_SF_vs_FBS.RData')


# get row_data:
row_data=count_wide%>%dplyr::select(Name,Description)
setnames(row_data,c('gene_id','gene_name'))


# format data: 
resOrdered=format_result(res,row_data)


# output rank for GSEA: 
output=resOrdered[,.(gene_name,stat)]
output=output[!is.na(stat),]
write.table(output,file='../processed_data/160801/hcasmc_fbs_vs_sf.rnk',quote=F,row.names=F,col.names=F,sep='\t')


# plot RPS16 (ENSG00000105193)
gene=data.frame(count=unlist(count_wide[Description=='RPS16',3:ncol(count_wide),with=F]))
gene=gene%>%mutate(condition=str_match(rownames(gene),'sf|fbs'))
p=ggplot(gene,aes(condition,count))+geom_boxplot()+ggtitle('RPS16')
save_plot('../figures/160801/RPS16.pdf',p)


# plot MRPL13
gene_name='MRPL13'
gene=data.frame(count=unlist(count_wide[Description==gene_name,3:ncol(count_wide),with=F]))
gene=gene%>%mutate(condition=str_match(rownames(gene),'sf|fbs'))
p=ggplot(gene,aes(condition,count))+geom_boxplot()+ggtitle(gene_name)
save_plot(sprintf('../figures/160801/%s.pdf',gene_name),p)
