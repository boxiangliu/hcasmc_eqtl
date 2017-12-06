gene_name=$1
snp=$2
out_prefix=$3
flank=${4:-500kb}
in_dir=${5:-/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/}
nikpay=${6:-/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/CARDIoGRAMplusC4D/cad.add.160614.website.metal.txt}
ukbb=${7:-/srv/persistent/bliu2/HCASMC_eQTL/data//gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz}

find_eQTL_file(){
	local gene_name=$1
	local in_dir=$2
	in_fn=$(ls $in_dir/chr*/*_$gene_name.pval.txt)
	
	local in_base=$(basename $in_fn)
	in_prefix=${in_base/.pval.txt/}
}
export -f find_eQTL_file


find_eQTL_file $gene_name $in_dir

echo INFO - $gene_name
echo INFO - $snp
echo INFO - $out_prefix
echo INFO - $in_fn
echo INFO - $in_prefix

mkdir -p $out_prefix/$in_prefix/


locuszoom \
--metal $in_fn \
--pvalcol pval --markercol rsid \
--refsnp $snp \
--flank $flank \
--source 1000G_March2012 --build hg19 --pop EUR \
title='eQTL' \
--no-date --prefix $out_prefix/$in_prefix/eQTL


locuszoom \
--metal $nikpay \
--pvalcol p_dgc --markercol markername \
--refsnp $snp \
--flank $flank \
--source 1000G_March2012 --build hg19 --pop EUR \
title='Nikpay' \
--no-date --prefix $out_prefix/$in_prefix//Nikpay


locuszoom \
--metal $ukbb  \
--pvalcol p-value_gc --markercol snptestid \
--refsnp $snp \
--flank $flank \
--source 1000G_March2012 --build hg19 --pop EUR \
title='UKBB' \
--no-date --prefix $out_prefix/$in_prefix/UKBB


