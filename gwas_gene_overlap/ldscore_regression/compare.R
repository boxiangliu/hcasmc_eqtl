library(data.table)

# Variables: 
baseline_ldscore_fn='/srv/persistent/bliu2/shared/ldscore/baseline/baseline.1.l2.ldscore.gz'
baseline_annot_fn='/srv/persistent/bliu2/shared/ldscore/baseline/baseline.1.annot.gz'

tissue_specific_ldscore_fn='../processed_data/gwas_gene_overlap/ldscore_regression/ldscore//Adipose - Subcutaneous.1.l2.ldscore.gz'
tissue_specific_annot_fn='../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_snp_annotation/Adipose - Subcutaneous.1.annot.gz'


# Read annot files: 
baseline_annot=fread(sprintf('zcat %s',baseline_annot_fn),select=1:3)
tissue_specific_annot=fread(sprintf('zcat "%s"',tissue_specific_annot_fn),select=1:3)
stopifnot(all(baseline_annot==tissue_specific_annot))


# Read LDscore files:
baseline_ldscore=fread(sprintf('zcat %s',baseline_ldscore_fn),select=1:3)
tissue_specific_ldscore=fread(sprintf('zcat "%s"',tissue_specific_ldscore_fn),select=1:3)
stopifnot(all(baseline_ldscore==tissue_specific_ldscore))


# merged=merge(baseline_ldscore,tissue_specific_ldscore,by=c('CHR','SNP','BP'),all=TRUE,sort=FALSE)

# zcat ../processed_data/gwas_gene_overlap/ldscore_regression/tissue_specific_snp_annotation/Adipose\ -\ Subcutaneous.1.annot.gz | wc -l 
# # 711594
# zcat $shared/ldscore/baseline/baseline.1.annot.gz | wc -l 
# # 711594
# zcat ../processed_data/gwas_gene_overlap/ldscore_regression/ldscore/Adipose\ -\ Subcutaneous.1.l2.ldscore.gz | wc -l
# # 98160
# zcat /srv/persistent/bliu2/shared/ldscore/baseline/baseline.1.l2.ldscore.gz | wc -l 
# # 98160