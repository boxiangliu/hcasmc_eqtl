library(data.table)
in_fn='../processed_data/gwas_atacseq_overlap/tmp/ld_set.tsv'
out_dir='../processed_data/gwas_atacseq_overlap/milos/gwas_proxy/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
ld_set=fread(in_fn)
fwrite(ld_set[loci_index==gwas_index],sprintf('%s/nikpay_gwas_proxy.tsv',out_dir),sep='\t')