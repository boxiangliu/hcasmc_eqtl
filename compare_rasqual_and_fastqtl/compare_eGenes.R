library(data.table)

ras_fn='../processed_data/rasqual/output_merged/treeQTL/eGenes.tsv'
fas_fn='../processed_data/eqtl/fastqtl/output/treeQTL/eGenes.tsv'

ras=fread(ras_fn)
fas=fread(fas_fn)

mean(fas$family%in%ras$family) # 0.8622754