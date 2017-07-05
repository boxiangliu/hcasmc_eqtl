library(data.table)
library(stringr)
library(gplots)

geno_fn='../processed_data/atacseq/asvcf/concordance/genotype.txt'
ase_fn='../processed_data/atacseq/asvcf/concordance/ase.txt'
fig_dir='../figures/atacseq/asvcf/concordance/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

geno=fread(geno_fn,select=2:9,header=FALSE)
geno=round(geno)


ase=fread(ase_fn,select=2:9,sep='\t',header=FALSE)
ref=apply(ase,2,FUN=function(x) {as.integer(str_split_fixed(x,',',2)[,1])})
alt=apply(ase,2,FUN=function(x) {as.integer(str_split_fixed(x,',',2)[,2])})
total=ref+alt
ratio=alt/total


concord=matrix(0,nrow=8,ncol=8)
rownames(concord)=colnames(concord)=c('1346','1508','1522','200212','2305','2356','2510','2989')
for (col1 in 1:8){

	row1=!is.na(ratio[,col1])
	
	for (col2 in 1:8){
	
		print(col1,col2)

		row2=unlist(geno[,col2,with=F])%in%c(0,2)
		row=row1&row2

		ratio_col=round(ratio[row,col1])
		geno_col=unlist(geno[row,col2,with=F])/2
		concord[col1,col2]=mean(geno_col==ratio_col)
	}
}

pdf(sprintf('%s/concordance.pdf',fig_dir))
heatmap.2(concord,Rowv='none',Colv='none',trace='none',dendrogram='none',xlab='ATACseq',ylab='WGS',margin=c(6,6))
dev.off()

