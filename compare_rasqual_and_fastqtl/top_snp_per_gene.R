library(data.table)
library(stringr)

out_dir='../processed_data/compare_rasqual_and_fastqtl/top_snp_per_gene/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


ras_fn='../processed_data/rasqual/output_merged/expressed_genes.pval.txt'

ras=fread(ras_fn)
ras[,c('hwe_chisq','impute_qual','delta','phi', 'overdispersion', 'sid', 'n_fsnp', 'n_rsnp', 'n_iter_null', 'n_iter_alt', 'tie', 'loglik_null', 'converge', 'r2_fsnp'):=NULL]
ras=ras[r2_rsnp>=0.8,]
ras[,min_pval:=min(pval),by='gene_id']
ras_top=ras[pval==min_pval,]
ras_top[,sid:=paste(chr,pos,ref,alt,'b37',sep='_')]
ras_top[,sid:=str_replace(sid,'chr','')]
fwrite(ras_sig,sprintf('%s/rasqual.txt',out_dir),sep='\t')
rm(ras)

fas_fn='../processed_data/eqtl/fastqtl/output/nominal/all.txt.gz'
fas=fread(sprintf('zcat %s',fas_fn))
setnames(fas,c('gene_id','sid','dist','pval','beta','se'))
fas[,min_pval:=min(pval),by='gene_id']
fas_top=fas[pval==min_pval,]
fwrite(fas_top,sprintf('%s/fastqtl.txt',out_dir),sep='\t')
rm(fas)