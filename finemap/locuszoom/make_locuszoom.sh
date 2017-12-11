gene_name=$1
snp=$2
out_prefix=$3
flank=${4:-500kb}
in_dir=${5:-/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/}
nikpay=${6:-/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/CARDIoGRAMplusC4D/cad.add.160614.website.metal.txt}
ukbb=${7:-/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz}
howson=${8:-/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/howson_2017/Howson-JMM_CHD_Mixed_2017.norm.in1kgp3.txt}
gtex_dir=${9:-/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_allpairs_FOR_QC_ONLY/}


find_eQTL_file(){
	local gene_name=$1
	local in_dir=$2
	in_fn=$(ls $in_dir/chr*/*_$gene_name.pval.txt)
	
	local in_base=$(basename $in_fn)
	in_prefix=${in_base/.pval.txt/}

	gene_id=${in_prefix/_$gene_name/}
}
export -f find_eQTL_file


find_eQTL_file $gene_name $in_dir

echo INFO - $gene_name
echo INFO - $snp
echo INFO - $out_prefix
echo INFO - $in_fn
echo INFO - $in_prefix

mkdir -p $out_prefix/$in_prefix/


# locuszoom \
# --metal $in_fn \
# --pvalcol pval --markercol rsid \
# --refsnp $snp \
# --flank $flank \
# --source 1000G_March2012 --build hg19 --pop EUR \
# title='eQTL' \
# --no-date --prefix $out_prefix/$in_prefix/eQTL


# locuszoom \
# --metal $nikpay \
# --pvalcol p_dgc --markercol markername \
# --refsnp $snp \
# --flank $flank \
# --source 1000G_March2012 --build hg19 --pop EUR \
# title='Nikpay' \
# --no-date --prefix $out_prefix/$in_prefix//Nikpay


# locuszoom \
# --metal $ukbb  \
# --pvalcol p-value_gc --markercol snptestid \
# --refsnp $snp \
# --flank $flank \
# --source 1000G_March2012 --build hg19 --pop EUR \
# title='UKBB' \
# --no-date --prefix $out_prefix/$in_prefix/UKBB


# locuszoom \
# --metal $howson  \
# --pvalcol p --markercol rsid \
# --refsnp $snp \
# --flank $flank \
# --source 1000G_March2012 --build hg19 --pop EUR \
# title='Howson' \
# --no-date --prefix $out_prefix/$in_prefix/Howson

for tissue in Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood; do
echo INFO - $tissue
bash finemap/locuszoom/find_rsid.sh \
$gene_id \
$gtex_dir/${tissue}_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz \
$out_prefix/$in_prefix/$tissue/fastQTL.txt


locuszoom \
--metal $out_prefix/$in_prefix/$tissue/fastQTL.txt  \
--pvalcol pval --markercol rsid \
--refsnp $snp \
--flank $flank \
--source 1000G_March2012 --build hg19 --pop EUR \
title=$tissue \
--no-date --prefix $out_prefix/$in_prefix/$tissue/$tissue
done
