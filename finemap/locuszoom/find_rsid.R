suppressWarnings(library(data.table))
suppressWarnings(library(stringr))

args=commandArgs(T)
in_fn=args[1]
out_fn=args[2]
message('INFO - input:', in_fn)
message('INFO - output:', out_fn)

message('INFO - reading eQTL file')
eqtl=fread(in_fn,col.names=c('gene_id','snp','dist','pval','beta','se'))
split_snp=str_split_fixed(eqtl$snp,'_',5)
eqtl$chr=split_snp[,1]
stopifnot(is.character(eqtl$chr))
eqtl$pos=as.integer(split_snp[,2])


start=min(eqtl$pos)
end=max(eqtl$pos)
chr=eqtl$chr[1]
message(sprintf('INFO - %s:%s-%s',chr,start,end))

tmp_fn=sprintf('%s.tmp',in_fn)
on.exit(unlink(tmp_fn))

vcf_fn=list.files('/mnt/lab_data/montgomery/shared/1KG/',pattern=sprintf('^ALL.chr%s\\..*.vcf.gz$',chr),full.names=TRUE)
cmd=sprintf('bcftools query -r %s:%s-%s -f "%%CHROM  %%POS  %%ID\n" %s > %s',chr,start,end,vcf_fn,tmp_fn)
message('INFO - command:', cmd)
system(cmd)

vcf=fread(tmp_fn,col.names=c('chr','pos','rsid'))
vcf$chr=as.character(vcf$chr)
stopifnot(is.integer(vcf$pos))
merged=merge(eqtl,vcf,by=c('chr','pos'))
fwrite(merged[,list(gene_id,rsid,dist,pval,beta,se)],out_fn,sep='\t')


