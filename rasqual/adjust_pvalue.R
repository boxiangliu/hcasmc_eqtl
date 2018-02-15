library(data.table)
library(stringr)
library(TreeQTL)

# command args: 
args=commandArgs(T)
in_file=args[1]
out_dir=args[2]
print(sprintf('INFO - input: %s',in_file))
print(sprintf('INFO - output dir: %s',out_dir))


# Function:
get_gene_map=function(genes){
	gencode=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')
	gencode[,geneid:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
	gene_map=unique(gencode[V3=='gene',.(geneid,chr=V1,left=V4,right=V5)])
	gene_map=gene_map[geneid%in%genes,]
	return(gene_map)
}

get_snp_map=function(snps){
	snp_map=data.table(snpid=unique(snps))
	snp_map[,chr:=str_split_fixed(snpid,'_',5)[,1]]
	snp_map[,chr:=paste0('chr',chr)]
	snp_map[,pos:=as.integer(str_split_fixed(snpid,'_',5)[,2])]
	return(snp_map)
}

treeQTL=function(meqtl,snp_map,gene_map,level1=0.05,level2=0.05,eSNP=TRUE,tmp_dir='./',cis_dist=1e6){

	fwrite(meqtl,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sep='\t')
	on.exit(unlink(sprintf('%s/mEQTL_out_cis.txt',tmp_dir)))
	on.exit(unlink(sprintf('%s/eAssoc_cis.txt',tmp_dir)))

	if (eSNP){
		print('INFO - level1 is eSNP')

		# Get number of genes nearby to each SNP:
		print('INFO - calculating number of tests per SNP...')
		n_tests_per_SNP=get_n_tests_per_SNP(snp_map,gene_map,nearby=TRUE,dist=cis_dist)

		# Get list of eSNPs:
		print('INFO - getting eSNPs...')
		eSNPs=get_eSNPs(n_tests_per_SNP, sprintf('%s/mEQTL_out_cis.txt',tmp_dir),level1=level1,level2=level2)
		setDT(eSNPs)

		# Generate txt output file with full set of eAssociations:
		print('INFO - getting eAssociations...')
		eAssoc=get_eAssociations(eSNPs,n_tests_per_SNP,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sprintf("%s/eAssoc_cis.txt",tmp_dir),by_snp=TRUE)
		setDT(eAssoc)

		return(list(eSNPs,eAssoc))

	} else {
		print('INFO - level1 is eGene')

		# Get number of SNPs nearby to each gene:
		print('INFO - calculating number of tests per gene...')
		n_tests_per_gene=get_n_tests_per_gene(snp_map, gene_map, nearby = TRUE, dist = cis_dist)


		# Get list of eGenes:
		print('INFO - getting eGenes...')
		eGenes=get_eGenes(n_tests_per_gene, sprintf('%s/mEQTL_out_cis.txt',tmp_dir),level1=level1,level2=level2)
		setDT(eGenes)

		# Generate txt output file with full set of eAssociations:
		print('INFO - getting eAssociations...')
		eAssoc=get_eAssociations(eGenes, n_tests_per_gene, sprintf('%s/mEQTL_out_cis.txt',tmp_dir), sprintf("%s/eAssoc_cis.txt",tmp_dir),by_snp=FALSE)
		setDT(eAssoc)

		return(list(eGenes,eAssoc))
	}
}



# read input: 
meqtl=fread(in_file,header=T)
setorder(meqtl,`p-value`)

gene_map=get_gene_map(unique(meqtl$gene))
snp_map=get_snp_map(unique(meqtl$SNP))


for (level in c(0.01,0.001)){
	# Correct p-value on eSNPs level:
	meqtl=meqtl[`p-value`<level,]
	res=treeQTL(meqtl,snp_map,gene_map,level1=level,level2=level,eSNP=TRUE)
	eSNPs=setDT(res[[1]])
	eAssoc=setDT(res[[2]])
	setorder(eSNPs,fam_p)
	setorder(eAssoc,BBFDR)

	fwrite(eSNPs,sprintf('%s/eSNPs_level%s.tsv',out_dir,level),sep='\t')
	fwrite(eAssoc,sprintf('%s/eAssoc_eSNPs_level%s.tsv',out_dir,level),sep='\t')


	# Correct p-value on eSNPs level:
	res=treeQTL(meqtl,snp_map,gene_map,level1=level,level2=level,eSNP=FALSE)
	eGenes=setDT(res[[1]])
	eAssoc=setDT(res[[2]])
	setorder(eGenes,fam_p)
	setorder(eAssoc,BBFDR)

	fwrite(eGenes,sprintf('%s/eGenes_level%s.tsv',out_dir,level),sep='\t')
	fwrite(eAssoc,sprintf('%s/eAssoc_eGenes_level%s.tsv',out_dir,level),sep='\t')
}

