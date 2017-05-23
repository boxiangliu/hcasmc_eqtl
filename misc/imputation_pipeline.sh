#!/bin/bash
# bosh liu
# 2016/04/10

input=$1

SHARED=/srv/persistent/bliu2/shared/
declare -A sizes=(["chr1"]=249250621 \
["chr2"]=243199373 \
["chr3"]=198022430 \
["chr4"]=191154276 \
["chr5"]=180915260 \
["chr6"]=171115067 \
["chr7"]=159138663 \
["chrX"]=155270560 \
["chr8"]=146364022 \
["chr9"]=141213431 \
["chr10"]=135534747 \
["chr11"]=135006516 \
["chr12"]=133851895 \
["chr13"]=115169878 \
["chr14"]=107349540 \
["chr15"]=102531392 \
["chr16"]=90354753 \
["chr17"]=81195210 \
["chr18"]=78077248 \
["chr20"]=63025520 \
["chrY"]=59373566 \
["chr19"]=59128983 \
["chr22"]=51304566 \
["chr21"]=48129895)

for chr in $(seq 1 22); do
	# subset to chromosome and generate bed files 
	bed=${input/.vcf.gzary/}
	bed=$bed.chr$chr.bed
	# echo "[ subset to chrosome ]"
	# echo "writing to $bed"
	# bash subset_to_chromosome.sh $input $bed $chr

	# check strand alignment:
	reference_hap=$SHARED/haplotype_reference/1000G_phase3/1000GP_Phase3/1000GP_Phase3_chr$chr.hap.gz
	reference_legend=$SHARED/haplotype_reference/1000G_phase3/1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz
	reference_sample=$SHARED/haplotype_reference/1000G_phase3/1000GP_Phase3/1000GP_Phase3.sample
	strand_alignment_output=${bed/.bed/.snp.strand}
	# echo "[ check strand alignment ]"
	# echo "using input:" $bed
	# echo "using reference panel:" $reference_hap $reference_legend $reference_sample
	# echo "writing strand alignment to:" $strand_alignment_output
	# bash check_strand_alignment.sh $bed $reference_hap $reference_legend $reference_sample $strand_alignment_output

	# plot percent missing from 1000G reference panel vs freq in main panel:
	freq=${bed/.bed/.frq}
	figure=../figures/plot_percentage_in_1000G.chr$chr.pdf
	echo "[ plot percentage missing from 1000G ]"
	echo "allele frequency file:" $freq
	echo "saving figure to:" $figure
	Rscript plot_percentage_in_1000G.R $freq $strand_alignment_output $figure

	# prephasing: 
	genetic_map=$SHARED/haplotype_reference/1000G_phase3/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt
	phased=${bed/.bed/.phased.haps}
	# echo "[ prephasing ]"
	# echo "excluding SNP in:" $strand_alignment_output.exclude
	# echo "using reference panel:" $reference_hap $reference_legend $reference_sample
	# echo "using genetic map:" $genetic_map
	# echo "writing phased haplotypes to:" $phased
	# bash prephasing.sh $bed $genetic_map $reference_hap $reference_legend $reference_sample $phased $strand_alignment_output.exclude

	# imputation: 
	imputed=${phased/.haps/.imputed}
	size=${sizes[chr$chr]}
	echo "[ imputation ]"
	echo "chr$chr has size $size"
	echo "using reference panel:" $reference_hap $reference_legend $reference_sample
	echo "using genetic map:" $genetic_map
	echo "writing imputed genotypes to:" $imputed
	bash impute.sh $phased $reference_hap $reference_legend $genetic_map $imputed $size
done 

touch .imputation_pipeline.sh.done