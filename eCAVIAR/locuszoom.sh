#!/bin/bash
# command args
gene=$1
# gene=ENSG00000236838.2
refsnp=$2
# refsnp=rs28524880
out_dir=$3
# out_dir=$figures/eCAVIAR/
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts
data=/srv/persistent/bliu2/HCASMC_eQTL/data

# subset to eQTL file to eGenes selected above: 
grep $gene ../processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt | sort -k1,2 -V > $out_dir/tmp.$gene

# add rs ID to each snp: 
cat $out_dir/tmp.$gene | awk 'BEGIN {FS="\t|_";OFS="\t"} {print $2,$3}' > $out_dir/tmp.${gene}.chr_pos
cat $out_dir/tmp.${gene}.chr_pos | python $scripts/160811/subset_dbsnp.py > $out_dir/tmp.${gene}.dbsnp
Rscript $scripts/160811/add_rsid.R -dbsnp=$out_dir/tmp.${gene}.dbsnp -eqtl=$out_dir/tmp.$gene -out=$out_dir/tmp.$gene.rsid


# make locuszoom plot for ENSG00000197208.5 at rs273909 (SLC22A4-SLC22A5):
cat $out_dir/tmp.$gene.rsid | awk 'BEGIN{OFS="\t"; print "markername","pval"} {print $7,$4}'> $out_dir/$gene.metal
locuszoom --metal $out_dir/$gene.metal --pvalcol pval --markercol markername --refsnp $refsnp --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL ($gene,$refsnp)" --prefix $out_dir/eqtl_${gene}


# cast CAD GWAS into METAL format: 
# cat $data/gwas/cad.add.160614.website.txt | awk 'BEGIN {OFS="\t"} {print $1, $11}' > $data/gwas/cad.add.160614.website.metal.txt

# make locuszoom plot for GWAS hit rs273909 (SLC22A4-SLC22A5):
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp $refsnp --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="GWAS ($gene,$refsnp)" --prefix $out_dir/gwas_${gene}
