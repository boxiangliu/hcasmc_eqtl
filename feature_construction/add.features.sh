#!/bin/bash

# author: Emily Tsang
# modified by Boxiang Liu to add features specifically for HCASMC eQTL

# INPUT
# * input bed file, see format below
# * 1000 genomes allele frequency file
#   (bed file with 5 additional columns with the AF for each super population [EAS, AMR, AFR, EUR, SAS])
# * output directory

#################################
# input bed file has
# 1* chromosome (with "chr")
# 2* Pos0
# 3* Pos1
# 4* Variant ID
# 5* Ref allele
# 6* Alt allele
# 7* MAF
# 8* CADD raw score
# 9* CADD phred score
#
# add the following annotations
# 10*  number of TFBS overlaps (Pouya's Roadmap motifs)
# 11* CpG (this and all features through 25 are from CADD)
# 12* priPhCons
# 13* mamPhCons
# 14* verPhCons
# 15* priPhyloP
# 16* mamPhyloP
# 17* verPhyloP
# 18* GerpN
# 19* GerpS
# 20* GerpRS
# 21* GerpRSpval
# 22* fitCons
# 23* TFBS
# 24* TFBSPeaks
# 25* TFBSPeaksMax
# 26* ATACseq peaks for sample 2305
# 27-31* 1KG allele frequency from 5 super populations (5 columns: EAS, AMR, AFR, EUR, SAS)
##################################

# REQUIRES
# * /mnt/lab_data/montgomery/shared/conservation/sorted.consolidated.annotations.bed.gz (15 cols, with chr)
# * /srv/scratch/restricted/goats/features/variantBeds/1KG/indels.1kg.AF.bed.gz (8 cols, with chr)
# * /srv/scratch/restricted/goats/features/variantBeds/1KG/SNPs.1kg.AF.bed.gz (ditto)
# * /srv/scratch/restricted/goats/features/annotations/TFBS_Pouya/*merged.bed.gz (3 cols, with chr)
# * /users/joed3/sardinia.eqtl.enrichments/data/er/e062.pbmcs.chromhmm.bed.gz 


# use development bedtools zach installed so that the sorted intersects are not an issue
export PATH=/users/zappala/software/BEDTools/bedtools2-dev/bedtools2/bin:$PATH

######################################

# Define the addFeatures function to run across all chromosomes

# PATHS TO SET #######################
annodir=/srv/scratch/restricted/goats/features/annotations # maybe change this for future version
tfdir=${annodir}/TFBS_Pouya
atac=/srv/scratch/bliu2/2305_ppr.IDR0.1.filt.narrowPeak.gz
consolidated=/mnt/lab_data/montgomery/shared/conservation/consolidated/sorted.consolidated.annotations.bed.gz
kgfile=/srv/scratch/restricted/goats/features/variantBeds/1KG/SNPs.1kg.AF.bed.gz

# Input and output filenames
f=$1
output_filename=$2
# Examples: 
# f=../data/variantBeds/chr1.withCADD
# output_filename=../data/variantBeds/chr1.allFeastures.bed.gz


# Add TF annotations
zcat $f.bed.gz | intersectBed -sorted -wa -c -a stdin -b ${tfdir}/*.merged.bed.gz > $f.withTF.bed

# Add consolidated conservation annotation
intersectBed -sorted -wa -wb -loj -a $f.withTF.bed -b ${consolidated} > $f.withTF.withCons.bed

# Add ER annotations
intersectBed -sorted -wa -wb -loj -a $f.withTF.withCons.bed -b ${atac} > $f.withTF.withCons.withATAC.bed

# Add 1KG MAFs
sed 's/\t\t/\t/g' $f.withTF.withCons.withATAC.bed | intersectBed -sorted -wa -wb -loj -a stdin -b ${kgfile} > $f.withTF.withCons.withATAC.with1KG.bed

# Process final output to simplified BED format
sed 's/\t\t/\t/g' $f.withTF.withCons.withATAC.with1KG.bed | awk 'BEGIN{tfEnd=10; consEnd=28; erEnd=32}{
	for(i=1; i<=tfEnd; i++){ printf("%s\t",$i)}; 
	for(i=tfEnd+4; i<=consEnd; i++){ printf("%s\t",$i)}; 
	for(i=consEnd+4; i<=erEnd; i++){ printf("%s\t",$i)}; 
	for(i=erEnd+4; i<NF; i++){ printf("%s\t",$i)};
	print $NF}' | sed 's/\t\./\tNA/g' | \
	gzip -c > $output_filename

# Remove intermediate files
rm $f.with*
