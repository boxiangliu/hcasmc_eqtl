# combine 2015 CAD and Metabochip variants: 
cat ../data/gwas/gwas_loci.cad.all.genomewide_fdr_merged.txt \
<(grep -v markername ../data/gwas/metabochip_novel_lead_variants.txt) > \
../data/gwas/gwas.cad.mi.recessive.metabochip.txt


# Select genes close to CAD GWAS hits (from the 2015 CAD and 2016 Metabochip variants): 
Rscript rasqual/select_genes.R 


# Run RASQUAL for all selected genes: 
grep -f ../processed_data/rasqual/prioritized_genes.txt -n ../processed_data/rasqual/input/rasqual.input.autosome.txt | awk '{print $1}' | awk 'BEGIN{FS=":"}{print $1}' > ../processed_data/rasqual/prioritized_genes.number.txt
parallel bash rasqual/rasqual.sh ../processed_data/rasqual/input/rasqual.input.autosome.txt {} ../processed_data/rasqual/Y.tidy.bin ../processed_data/rasqual/K.bin ../processed_data/rasqual/X.bin ../processed_data/genotype/phasing_with_1kg/phased_and_imputed_gprobs/phased_and_imputed.all.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz :::: ../processed_data/rasqual/prioritized_genes.number.txt 


# Run RASQUAL for other genes: 
parallel -j10 bash rasqual/rasqual.sh ../processed_data/rasqual/input/rasqual.input.autosome.txt {} ../processed_data/rasqual/Y.tidy.bin ../processed_data/rasqual/K.bin ../processed_data/rasqual/X.bin ../processed_data/genotype/phasing_with_1kg/phased_and_imputed_gprobs/phased_and_imputed.all.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz :::: ../processed_data/rasqual/unprioritized_genes.number.txt 

