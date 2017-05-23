library(data.table)
library(stringr)
library(cowplot)


# Variables: 
in_fn='../processed_data/gwas_atacseq_overlap/gregor/merge_peaks_adult_filt/HCASMC.merged.bed'
fig_dir='../figures/eqtl_and_atacseq/gregor/overlap_enrichment/'
out_dir='../processed_data/eqtl_and_atacseq/gregor/overlap_enrichment/'
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
		s=sum(overlap[gwas_index==snpID,loci_overlap])
		p4_=sum(S_sum>=s)/length(S_sum)
	} else {
		p4_=p4(u_hat,n,p)
	}
	return(list(u_hat=u_hat,p4=p4_))
}

# Read GWAS and matched background SNPs (and LD SNPs):
ld_set=fread('../processed_data/eqtl_and_atacseq/tmp/background_set.top1000eqtl.tsv')
ld_set[,c('start','end'):=list(pos,pos)]
ld_set[,chr:=paste0('chr',chr)]
setkey(ld_set,chr,start,end)


# Calculate enrichment statistics for only adult tissue/cell line (and remove AUDIT ERROR samples): 
# Read HCASMC ATACseq:
dhs=fread(in_fn)
dhs[,c('start','end'):=list(start-500,end+500),]
setkey(dhs,chr,start,end)


# Overlap: 
overlap=unique(foverlaps(ld_set,dhs[,list(chr,start,end)]))
overlap[,c('i.start','i.end'):=NULL]


overlap[,loci_overlap:=!is.na(start)]
overlap[,c('start','end'):=NULL]
overlap=unique(overlap)
overlap[,p:=mean(loci_overlap),by='eqtl']


# Calculate enrichment p-value:
p=overlap[,list(p=unique(p)),by='eqtl']
n=rep(1,length(p$p))
s=sum(overlap[snpID==eqtl,loci_overlap])
pval=p_non_identical_binom(n,p$p,s)$p4
# 3.03409e-05