# Select LD variants
# Boxiang Liu (bliu2@stanford.edu)
# 2017-11-15
library(foreach)
library(doMC)
registerDoMC(10)
library(data.table)

select_LD_variants=function(index_set,tmp_dir,vcf_dir){
	# Calculate LD r2:
	foreach(c=1:22)%dopar%{
		ld_out_prefix=sprintf('%s/ld0.7.chr%s',tmp_dir,c)

		snp_list_fn=sprintf('%s/chr%s.snps',tmp_dir,c)
		message('INFO - writing SNP list to ',snp_list_fn)
		write.table(data.frame(unique(index_set[chr==c,snpID])),file=snp_list_fn,col.names=FALSE,row.names=FALSE,quote=FALSE)

		message('INFO - calculating LD for chr',c)
		command=sprintf('plink --vcf %s/chr%s.vcf.gz --keep-allele-order --r2 --ld-snp-list %s --ld-window-kb 1000 --ld-window-r2 0.7 --out %s',vcf_dir,c,snp_list_fn,ld_out_prefix)
		print(command)
		system(command)
	}
	command=sprintf('cat %s/ld0.7.chr*.ld | grep -v "CHR_A" > %s/ld0.7.all_chr.ld',tmp_dir,tmp_dir)
	system(command)


	# Get Create LD snp set:
	message('INFO - creating LD snp set...')
	ld=fread(sprintf('%s/ld0.7.all_chr.ld',tmp_dir),col.names=c('CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R2'))

	ld_set=foreach(i=1:nrow(index_set),.combine='rbind')%dopar%{
		background_variant=index_set[i,snpID]
		foreground_variant=index_set[i,foreground_variant]

		ld_snp=data.table(snpID=ld[SNP_A==background_variant,SNP_B],r2=ld[SNP_A==background_variant,R2])
		ld_snp$background_variant=background_variant
		ld_snp$foreground_variant=foreground_variant
		ld_snp[,ld_proxy:=snpID!=background_variant]
		return(ld_snp)
	}

	ld_set[,chr:=paste0('chr',str_split_fixed(snpID,':',2)[,1])]
	ld_set[,pos:=as.integer(str_split_fixed(snpID,':',2)[,2])]
	return(ld_set)
}

