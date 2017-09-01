library(data.table)
library(stringr)
library(cowplot)
library(gridExtra)

# Variables: 
in_dir='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks/'
in_dir_filt='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_filt/'
in_dir_adult='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult/'
in_dir_adult_filt='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult_filt/'
in_dir_2007_2012='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_2007_2012/'
in_dir_2007_2012_noCancer='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_2007_2012_noCancer/'
fig_dir='../figures/gwas_atacseq_overlap/gregor/overlap_enrichment/'
out_dir='../processed_data/gwas_atacseq_overlap/gregor/overlap_enrichment/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Functions: 
p_non_identical_binom=function(n,p,s){
	# using saddle point approximation from Te Grotenhuis 2013
	# solve saddle point: 
	q=function(u,p){
		return(p*exp(u)/(1-p+p*exp(u)))
	}
	K=function(u,n,p){
		return(sum(n*log(1-p+p*exp(u))))
	}
	Kp=function(u,n,p){
		return(sum(n*q(u,p)))
	}
	Kpp=function(u,n,p){
		q_=q(u,p)
		return(sum(n*q_*(1-q_)))
	}
	Kppp=function(u,n,p){
		q_=q(u,p)
		return(sum(n*q_*(1-q_)*(1-2*q_)))
	}
	Kpppp=function(u,n,p){
		q_=q(u,p)
		return(sum(n*q_*(1-q_)*(1-6*q_*(1-q_))))
	}
	saddlepoint=function(u,n,p,s){
		return(Kp(u,n,p)-s)
	}
	w=function(u,n,p){
		K_=K(u,n,p)
		Kp_=Kp(u,n,p)
		return(sign(u)*sqrt(2*u*Kp_-2*K_))
	}
	u1=function(u,n,p){
		Kpp_=Kpp(u,n,p)
		return((1-exp(-u))*sqrt(Kpp_))
	}
	p3=function(u,n,p){
		w_=w(u,n,p)
		u1_=u1(u,n,p)
		if (isTRUE(all.equal(w_,0)) & isTRUE(all.equal(u1_,0))) {
			return(1/2-(2*pi)^(-1/2)*(1/6)*Kppp(0,n,p)*K(0,n,p)^(-3/2)-(1/2)*Kpp(0,n,p)^(-1/2))
		} else {
			return(1-pnorm(w_)-dnorm(w_)*(1/w_-1/u1_))
		}
	}

	u2=function(u,n,p){
		return(u*sqrt(Kpp(u,n,p)))
	}
	k3=function(u,n,p){
		return(Kppp(u,n,p)*Kpp(u,n,p)^(-3/2))
	}
	k4=function(u,n,p){
		return(Kpppp(u,n,p)*Kpp(u,n,p)^(-2))
	}
	p4=function(u,n,p){
		u2_=u2(u,n,p)
		k3_=k3(u,n,p)
		k4_=k4(u,n,p)
		w_=w(u,n,p)
		p3_=p3(u,n,p)
		return(p3_-dnorm(w_)*((1/u2_)*((1/8)*k4_-(5/24)*k3_^2)-1/(u2_^3)-k3_/(2*u2_^2)+1/w_^3))
	}

	u_hat=tryCatch(
		{uniroot(saddlepoint,lower=-100,upper=100,tol = 0.0001,n=n,p=p,s=s)$root},
		error=function(w){
			return(NA)
		})
	if (is.na(u_hat)){
		n_exp=100000
		n_gwas_index=length(p)
		set.seed(42)
		mat=matrix(rbinom(n_exp*n_gwas_index,1,p),nrow=n_gwas_index,ncol=n_exp)
		S_sum=colSums(mat)
		# s=sum(overlap[gwas_index==snpID,loci_overlap])
		p4_=sum(S_sum>=s)/length(S_sum)
	} else {
		p4_=p4(u_hat,n,p)
	}
	return(list(u_hat=u_hat,p4=p4_))
}

calc_enrichment_stat=function(peak_dir,ld_set){
	fn=list.files(peak_dir,pattern='bed')
	pval=c()
	tissue=c()
	for (f in fn){
		sample=str_replace(f,'.merged.bed','')
		tissue=c(tissue,sample)
		print(sprintf('INFO - %s',sample))
		dhs=fread(sprintf('%s/%s',peak_dir,f),col.names=c('chr','start','end'))
		setkey(dhs,chr,start,end)


		# Overlap: 
		overlap=unique(foverlaps(ld_set,dhs[,list(chr,start,end)]))
		overlap[,c('i.start','i.end'):=NULL]


		overlap[,snp_overlap:=!is.na(start)]
		overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
		overlap[,c('start','end'):=NULL]
		overlap=unique(overlap)
		stopifnot(nrow(overlap)==nrow(ld_set))
		stopifnot(overlap$snpID==ld_set$snpID)


		overlap=overlap[ld_proxy==FALSE,]
		overlap[,p:=mean(loci_overlap),by='gwas_index']


		# Calculate enrichment p-value:
		p=overlap[,list(p=unique(p)),by='gwas_index']
		n=rep(1,length(p$p))
		s=sum(overlap[snpID==gwas_index,loci_overlap])
		print(nrow(overlap))
		p_non_identical_binom(n,p$p,s)
		print('1')
		pval=c(pval,p_non_identical_binom(n,p$p,s)$p4)
		print('2')
	}
	pval=data.table(tissue,pval)
	setorder(pval,pval)
	return(pval)
}


# Read GWAS and matched background SNPs (and LD SNPs):
ld_set=fread('../processed_data/gwas_atacseq_overlap/tmp/ld_set.tsv')
ld_set[,c('start','end'):=list(pos,pos)]
setkey(ld_set,chr,start,end)


# Calculate enrichment statistics for all tissue/cell line: 
pval=calc_enrichment_stat(in_dir,ld_set)
pdf(sprintf('%s/gregor_pval_all_life_stages.pdf',fig_dir));grid.table(head(pval,20));dev.off()
fwrite(pval,sprintf('%s/gregor_pval_all_life_stages.tsv',out_dir),sep='\t')


# Plot correlation between p-value and life stage:
metadata=fread('../data/encode/dnase_seq/metadata.tsv')
metadata=metadata[`Biosample type`%in%c('tissue','primary cell')]

tmp=metadata[,list(adult_pct=mean(`Biosample life stage`=='adult'),fetus_pct=mean(`Biosample life stage`=='fetal')),by='Biosample term name']
tmp=rbind(tmp,data.table(`Biosample term name`='HCASMC',adult_pct=1,fetus_pct=0))

merged=merge(pval,tmp,by.x='tissue',by.y='Biosample term name',sort=F)
p1=ggplot(merged[(adult_pct==1)|(fetus_pct==1),],aes(pval,fill=ifelse(adult_pct==1,'adult','fetal')))+geom_histogram(position=position_dodge(width=0.15),binwidth=0.2)+scale_fill_discrete(name='Life stage')
save_plot(sprintf('%s/adult_vs_fetal.pdf',fig_dir),p1,base_height=6)


# Calculate enrichment statistics for all tissue/cell line (but without audit error):
pval=calc_enrichment_stat(in_dir_filt,ld_set)
pdf(sprintf('%s/gregor_pval_all_life_stages_filt.pdf',fig_dir));grid.table(head(pval,20));dev.off()
fwrite(pval,sprintf('%s/gregor_pval_all_life_stages_filt.tsv',out_dir),sep='\t')


# Calculate enrichment statistics for only adult tissue/cell line: 
fn=list.files(in_dir_adult,pattern='bed')
pval=c()
tissue=c()
for (f in fn){
	sample=str_replace(f,'.merged.bed','')
	tissue=c(tissue,sample)
	print(sprintf('INFO - %s',sample))
	dhs=fread(sprintf('%s/%s',in_dir_adult,f))
	setkey(dhs,chr,start,end)


	# Overlap: 
	overlap=unique(foverlaps(ld_set,dhs[,list(chr,start,end)]))
	overlap[,c('i.start','i.end'):=NULL]


	overlap[,snp_overlap:=!is.na(start)]
	overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
	overlap[,c('start','end'):=NULL]
	overlap=unique(overlap)
	stopifnot(nrow(overlap)==nrow(ld_set))


	overlap=overlap[ld_proxy==FALSE,]
	overlap[,p:=mean(loci_overlap),by='gwas_index']


	# Calculate enrichment p-value:
	p=overlap[,list(p=unique(p)),by='gwas_index']
	n=rep(1,length(p$p))
	s=sum(overlap[snpID==gwas_index,loci_overlap])
	pval=c(pval,p_non_identical_binom(n,p$p,s)$p4)
}
pval=data.table(tissue,pval)
setorder(pval,pval)
pdf(sprintf('%s/gregor_pval_adult.pdf',fig_dir));grid.table(head(pval,20));dev.off()
fwrite(pval,sprintf('%s/gregor_pval_adult.tsv',out_dir),sep='\t')

# Calculate enrichment statistics for only adult tissue/cell line (and remove AUDIT ERROR samples): 
fn=list.files(in_dir_adult_filt,pattern='bed')
pval=c()
tissue=c()
for (f in fn){
	sample=str_replace(f,'.merged.bed','')
	tissue=c(tissue,sample)
	print(sprintf('INFO - %s',sample))
	dhs=fread(sprintf('%s/%s',in_dir_adult_filt,f))
	setkey(dhs,chr,start,end)


	# Overlap: 
	overlap=unique(foverlaps(ld_set,dhs[,list(chr,start,end)]))
	overlap[,c('i.start','i.end'):=NULL]


	overlap[,snp_overlap:=!is.na(start)]
	overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
	overlap[,c('start','end'):=NULL]
	overlap=unique(overlap)
	stopifnot(nrow(overlap)==nrow(ld_set))


	overlap=overlap[ld_proxy==FALSE,]
	overlap[,p:=mean(loci_overlap),by='gwas_index']


	# Calculate enrichment p-value:
	p=overlap[,list(p=unique(p)),by='gwas_index']
	n=rep(1,length(p$p))
	s=sum(overlap[snpID==gwas_index,loci_overlap])
	pval=c(pval,p_non_identical_binom(n,p$p,s)$p4)
}
pval=data.table(tissue,pval)
setorder(pval,pval)
pdf(sprintf('%s/gregor_pval_adult_filt.pdf',fig_dir));grid.table(head(pval,20));dev.off()
fwrite(pval,sprintf('%s/gregor_pval_adult_filt.tsv',out_dir),sep='\t')
pval[,tissue:=factor(tissue,levels=tissue)]
p1=ggplot(pval,aes(tissue,-log10(pval)))+geom_point()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
save_plot(sprintf('%s/gregor_pval_adult_filt.scatter.pdf',fig_dir),p1,base_height=6)

# Calculate enrichment statistics for uniformly processed cell lines from 2007-2012
fn=list.files(in_dir_2007_2012,pattern='bed')
pval=c()
tissue=c()
for (f in fn){
	sample=str_replace(f,'.merged.bed','')
	tissue=c(tissue,sample)
	print(sprintf('INFO - %s',sample))
	dhs=fread(sprintf('%s/%s',in_dir_2007_2012,f))
	setkey(dhs,chr,start,end)


	# Overlap: 
	overlap=unique(foverlaps(ld_set,dhs[,list(chr,start,end)]))
	overlap[,c('i.start','i.end'):=NULL]


	overlap[,snp_overlap:=!is.na(start)]
	overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
	overlap[,c('start','end'):=NULL]
	overlap=unique(overlap)
	stopifnot(nrow(overlap)==nrow(ld_set))


	overlap=overlap[ld_proxy==FALSE,]
	overlap[,p:=mean(loci_overlap),by='gwas_index']


	# Calculate enrichment p-value:
	p=overlap[,list(p=unique(p)),by='gwas_index']
	n=rep(1,length(p$p))
	s=sum(overlap[snpID==gwas_index,loci_overlap])
	pval=c(pval,p_non_identical_binom(n,p$p,s)$p4)
}
pval=data.table(tissue,pval)
setorder(pval,pval)
pdf(sprintf('%s/gregor_pval_2007_2012.pdf',fig_dir));grid.table(head(pval,20));dev.off()
fwrite(pval,sprintf('%s/gregor_pval_2007_2012.tsv',out_dir),sep='\t')


# Calculate enrichment statistics for 
# uniformly processed cell lines from 2007-2012 (no cancer lines)
pval=calc_enrichment_stat(in_dir_2007_2012_noCancer,ld_set)
pdf(sprintf('%s/gregor_pval_2007_2012_noCancer.pdf',fig_dir));grid.table(head(pval,20));dev.off()
fwrite(pval,sprintf('%s/gregor_pval_2007_2012_noCancer.tsv',out_dir),sep='\t')


# Count the number of SNPs falling into each tissue/cell line:
fn=list.files(in_dir_adult_filt,pattern='bed')
n_overlap=c()
tissue=c()
gwas_set=ld_set[snpID==gwas_index]
setkey(gwas_set,chr,start,end)
for (f in fn){
	sample=str_replace(f,'.merged.bed','')
	tissue=c(tissue,sample)
	print(sprintf('INFO - %s',sample))
	dhs=fread(sprintf('%s/%s',in_dir_adult_filt,f))
	setkey(dhs,chr,start,end)


	# Overlap: 
	overlap=unique(foverlaps(gwas_set,dhs[,list(chr,start,end)]))
	overlap[,c('i.start','i.end'):=NULL]


	overlap[,snp_overlap:=!is.na(start)]
	overlap[,loci_overlap:=any(snp_overlap),by='loci_index']
	overlap[,c('start','end'):=NULL]
	overlap=unique(overlap)
	stopifnot(nrow(overlap)==nrow(gwas_set))


	overlap=overlap[ld_proxy==FALSE,]
	n_overlap=c(n_overlap,unlist(overlap[,list(n=sum(loci_overlap))]))
}

n_overlap=data.table(tissue,n_overlap)
setorder(n_overlap,-n_overlap)
pdf(sprintf('%s/gregor_num_overlap_adult_filt.pdf',fig_dir));grid.table(head(n_overlap,20));dev.off()
fwrite(n_overlap,sprintf('%s/gregor_num_overlap_adult_filt.tsv',out_dir),sep='\t')
