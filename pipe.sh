###################################
# Pipeline for HCASCM project     #
# Author: Boxiang Liu 			  #
# Contact: jollier.liu@gmail.com  #
###################################


HG19 = /srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
GATK =  /usr/bin/GenomeAnalysisTK.jar
GRCh37 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa
SHARED = /srv/persistent/bliu2/shared/

#--------- Genotypes -------------#
# Setup: 
mkdir genotype ../figures/genotype


# Prepare data: 
bash genotype/preprocessing/copy_WGS_HCASMC_from_valk.sh
python genotype/preprocessing/compress_gVCF.py /mnt/data/WGS_HCASMC sample_list.txt
python genotype/preprocessing/check_finished_sample.py /mnt/data/WGS_HCASMC sample_list.txt
python genotype/preprocessing/copy_bam_and_vcf_to_valk.sh /mnt/data/WGS_HCASMC/sample_list.txt
python genotype/preprocessing/concat_1000G_vcfs.py ../../shared/1000genomes/phase3v5a sample_list.txt ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf
bash genotype/preprocessing/fasta_dict.sh


# Joint genotyping and recalibration: 
mkdir ../data/joint2
bash genotype/joint_genotyping/genotype_gvcfs.sh
bash genotype/joint_genotyping/recalibrate_SNP.sh raw_variants.vcf
bash genotype/joint_genotyping/recalibrate_INDEL.sh recalibrated_snps_raw_indels.vcf
bash genotype/joint_genotyping/recalibrate_cleanup.sh 
Rscript genotype/joint_genotyping/plot_sensitivity.R # seems like 98% sensitivity (elbow) is a good cutoff. 
bash genotype/joint_genotyping/filter_variants.sh
bash genotype/joint_genotyping/biallelic_snp.sh
bash genotype/joint_genotyping/rename_chr.sh


# Quality control:
bash genotype/quality_control/vcf_stats.sh
grep "^SN" ../processed_data/genotype/quality_control/stats.txt | awk 'BEGIN {FS="\t"} {print $0}' > ../processed_data/genotype/quality_control/count.txt
Rscript genotype/quality_control/variant_count_by_type.R


# Phasing (and imputation) without reference panels: 
bash genotype/phasing_no_ref/download.sh
bash genotype/phasing_no_ref/beagle_no_ref.sh


# Post imputation quality control:
zgrep -v "^#" ../data/joint2/recalibrated_biallelic_SNP.beagle.vcf.gz | awk 'BEGIN {FS="\t|;|="; OFS="\t"; print "CHROM","POS","AR2","DR2","AF"} {print $1,$2,$9,$11,$13}' >  ../processed_data/genotype/phasing_no_ref/beagle_QC/recalibrated_biallelic_SNP.r2.tsv
bash genotype/phasing_no_ref/beagle_QC.R -input=../processed_data/genotype/phasing_no_ref/beagle_QC/recalibrated_biallelic_SNP.r2.tsv -figure_dir=../figures/160515_beagle_QC/
bash genotype/phasing_no_ref/update_sample_names.sh
Rscript genotype/quality_control/gen_sample_sheet_each_ethnicity.R # output Caucasian.txt, Asian.txt, AA.txt, Hispanic.txt
bash genotype/quality_control/run_detect_WGS_contamination.sh


# Post imputation filtering:
bash genotype/phasing_no_ref/filter_r2.sh
bash genotype/phasing_no_ref/filter_hwe.sh
rsync -vzh /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.hwe.vcf bosh@valkyr.stanford.edu:/home/diskstation/wgs/WGS_HCASMC_working_data_set/
Rscript genotype/phasing_no_ref/make_sample_list.R
bash genotype/phasing_no_ref/extract_dosage.sh


# Genotype PC:
bash genotype/genotype_pc/genotype_pc.sh


# Phasing (and imputation) with 1000 Genome Reference Panel: 
bash genotype/phasing_with_1kg/conform_gt.sh
bash genotype/phasing_with_1kg/phasing_with_1kg.sh


# Merge 15 array with 52 WGS data:
bash genotype/merge52and15/merge.ar2-0.4.sh 

#--------------- RNAseq ----------------# 
# Setup: 
mkdir ../data/rnaseq2


# Preprocessing:
bash rnaseq/preprocess/get_rnaseq_data.sh
bash rnaseq/preprocess/sort_and_index.sh
bash rnaseq/preprocess/quickcheck.sh
bash rnaseq/preprocess/RNAseQC.sh
bash rnaseq/preprocess/copy_rpkm.sh
bash rnaseq/preprocess/combine_rpkm.sh
Rscript rnaseq/preprocess/filter_rpkm.R 0.1 10 ../processed_data/rnaseq/preprocess/combine_rpkm/combined.rpkm ../processed_data/rnaseq/preprocess/combine_rpkm/combined.filter.rpkm


# Quality control:
Rscript rnaseq/quality_control/sample_correlation.R


#----------------- RNA-WGS match -------------#
bash rna_wgs_match/hg19toGRCh37.sh
bash rna_wgs_match/verifyBamID.sh 


#----------------- eQTL ------------------------
# Setup:
mkdir eqtl ../processed_data/eqtl ../figures/eqtl

# extract PEER factors:
bash eqtl/peer/peer_factor.sh
Rscript eqtl/peer/peer_factor_correlation.R


# Matrix eQTL:
bash eqtl/matrix_eqtl/combine_covariates.sh
Rscript eqtl/matrix_eqtl/covariates_correlation.R
bash eqtl/matrix_eqtl/matrix_eqtl.sh


# FastQTL:
bash eqtl/fastqtl/change_sid.sh
bash eqtl/fastqtl/fastqtl.nominal.sh
bash eqtl/fastqtl/adjust_pvalue.nominal.sh

#------------------ RASQUAL ---------------------
# Setup:  
mkdir rasqual ../figures/rasqual ../processed_data/rasqual


# Preprocessing:
bash rasqual/count_table.sh
python rasqual/calc_gcc.py /mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa /srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf exon > ../processed_data/rasqual/expression/gcc.exon.txt
bash rasqual/calc_offset.sh
bash rasqual/vcf2asvcf.wrapper.sh
bash rasqual/make_input.sh
Rscript $scripts/rasqual/combine_covariates.R --genotype_pc=../processed_data/genotype/genotype_pc/genotype_pcs.52samples.tsv --peer=../processed_data/eqtl/peer/factors.tsv --sample_info=../data/sample_info/sample_info.xlsx --output=../processed_data/rasqual/X.txt --gender_coding=numerical --num_geno_pc=4 --num_peer_factor=8
R --vanilla --quiet --args ../processed_data/rasqual/Y.tidy.txt ../processed_data/rasqual/K.txt ../processed_data/rasqual/X.txt < ../../tools/rasqual/R/txt2bin.R 


# Run RASQUAL:
bash rasqual/rasqual.wrapper.sh
bash rasqual/calc_pval.wrapper.sh
bash rasqual/merge_output.sh
bash rasqual/adjust_pvalue.sh



# concatenate all priortized genes: 
bash rasqual/rasqual_output2eCAVIAR_input.sh \
../processed_data/rasqual/prioritized_genes.output_file.txt \
../processed_data/rasqual/merged_output/prioritized_genes.allpairs.txt



# make locuszoom for TCF21 at rs2327429:
mkdir ../figures/rasqual/locuszoom/
cut -f2,26 ../processed_data/rasqual/output/ENSG00000118526.6_TCF21_w2000000.pval.txt > ../processed_data/rasqual/output/ENSG00000118526.6_TCF21_w2000000.metal
locuszoom --metal ../processed_data/rasqual/output/ENSG00000118526.6_TCF21_w2000000.metal --pvalcol pval --markercol rsid --refsnp rs2327429 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (TCF21)" --prefix ../figures/rasqual/locuszoom/TCF21_eQTL


# make locuszoom for TCF21 at rs12190287:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000118526.6_TCF21_w2000000.metal --pvalcol pval --markercol rsid --refsnp rs12190287 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (TCF21)" --prefix ../figures/rasqual/locuszoom/TCF21_eQTL


# make locuszoom for GWAS at rs12190287:
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs12190287 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="GWAS (rs12190287)" --prefix ../figures/rasqual/locuszoom/TCF21_gwas



#---------------- Compare RASQUAL and fastQTL ------------# 
# Setup: 
mkdir compare_rasqual_and_fastqtl ../processed_data/compare_rasqual_and_fastqtl

# Subset the gene, SNP, p-value, q-value, and effect size column of fastQTL and RASQUAL result: 
bash compare_rasqual_and_fastqtl/subset.sh

# Calculate the percentage of RASQUAL eQTL also discovered fastQTL: 
bash compare_rasqual_and_fastqtl/compare.R

# Replication in GTEx:
bash compare_rasqual_and_fastqtl/define_gtex_eqtl.sh
bash compare_rasqual_and_fastqtl/gtex_replication.R

# Overlap with ATACseq:
bash compare_rasqual_and_fastqtl/top_snp_per_gene.R
bash compare_rasqual_and_fastqtl/atacseq_overlap.R


#----- trans eQTL -------
# run matrix eQTL to map trans eQTL:
covariates=$processed_data/find_optimum_num_PEER_factors_matrixeqtl/covariates.matrixeqtl.pc3.peer8.tsv
Rscript $scripts/run_matrix_eQTL.R \
	$processed_data/dosage.tsv \
	$processed_data/genotype_loc.txt \
	$processed_data/combined.filter.norm.rpkm \
	$processed_data/gene_loc.txt \
	$covariates \
	$processed_data/cutoff0.05. \
	0 \
	0.05


# prepare data for Sherlock:
mkdir $processed_data/sherlock
awk 'BEGIN {OFS="\t"} {print $2,$1,$5}' $processed_data/trans.txt > $processed_data/sherlock/trans.sherlock.txt
awk 'BEGIN {OFS="\t"} {print $2,$1,$5}' $processed_data/trans.txt > $processed_data/sherlock/trans.sherlock.txt
zcat $processed_data/find_optimum_num_PEER_factors_matrixeqtl/pc3.peer8.2.cis.txt.gz | awk 'BEGIN {OFS="\t"} {print $2,$1,$5}' > $processed_data/sherlock/cis.sherlock.txt

# 16/06/03:
# setup: 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603
mkdir $scripts $processed_data $figures


# copy dASE RNAseq alignements from valk: 
bash $scripts/copy_rnaseq_from_valk.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments/sample_list.txt

# rename dASE RNAseq samples: 
bash $scripts/rename_rnaseq_samples.sh \
	$scripts/rename_rnaseq_samples.map.txt \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments/

# run RNAseQC: 
cp 160517_run_RNAseQC.sh 160603/run_RNAseQC.sh
bash $scripts/run_RNAseQC.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments


# combine rpkms from each sample to two files,
# one for fbs and one for sf samples 
cp 160527/combine_rpkm.sh $scripts
bash $scripts/combine_rpkm.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/rpkm \
	$processed_data/rpkm


# copy rpkm to the rpkm folder:
bash $scripts/copy_rpkm.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/rpkm


# prepare GTEx and hcasmc rpkm files (all genes): 
mkdir $processed_data/rpkm
bash $scripts/prepare_gtex_rpkm.sh $processed_data/rpkm
ln ../processed_data/160527/combined.rpkm $processed_data/rpkm/hcasmc.rpkm


# combine GTEx and HCASMC rpkm files: 
# subset to > 0.1 rpkm > 10 individuals: 
# log2(x+2) transform
Rscript $scripts/combine_and_filter_rpkm.R \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/sample_list.txt \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined \
	0.1 \
	10 

# hierarchical clustering:
Rscript $scripts/hclust.R \
	-rpkm=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm \
	-coldata=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col \
	-figure=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/hclust.pdf

# multidimensional scaling (2D):
Rscript $scripts/mds.R \
	-rpkm=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm \
	-coldata=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col \
	-tissue_names=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.txt \
	-figure=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/

# multidimensional scaling (3D).
# Since rgl is not installed on durga, the script needs to be run locally. 
# Rscript $scripts/mds.3D.R \
# 	-rpkm=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm \
# 	-coldata=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col \
# 	-tissue_names=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.txt \
# 	-figure=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/





# 16/06/14
# setup: 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160614
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160614
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160614
mkdir $scripts $processed_data $figures


# generate read counts: 
cp htseq_count.sh $scripts/
bash $scripts/htseq_count.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments \
	$scripts/htseq.eqtl.sample_list.txt \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/read_count

bash $scripts/htseq_count.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/alignments \
	$scripts/htseq.dase.sample_list.txt \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq_dase/read_count


# get read count from RNAseQC output: 
Rscript $scripts/rnaseqc_count.R \
	-input=$scripts/rnaseqc_count.sample_dirs.eqtl.txt

Rscript $scripts/rnaseqc_count.R \
	-input=$scripts/rnaseqc_count.sample_dirs.dase.txt


# merge rnaseqc output: 
Rscript $scripts/merge_rnaseqc_output.R \
	-input=$scripts/merge_rnaseqc_output.eqtl.txt \
	-output=$processed_data/rnaseqc.hcasmc_eqtl.reads.gct
Rscript $scripts/merge_rnaseqc_output.R \
	-input=$scripts/merge_rnaseqc_output.dase.txt \
	-output=$processed_data/rnaseqc.hcasmc_dase.reads.gct


 # differential expression against fibroblast: 
 Rscript $scripts/DESeq2.fibroblast.R


# 16/06/15: 
# setup: 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160615
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160615
mkdir $scripts $processed_data $figures


# download GWAS p-values: 
wget -P $processed_data/gwas http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip
wget -P $processed_data/gwas http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.recessive.Oct2015.pub.zip
wget -P $processed_data/gwas http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/mi.additive.Oct2015.pub.zip


# get top eQTLs from HCASMC:
Rscript $scripts/get_top_eqtls.hcasmc.R \
-input=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/find_optimum_num_PEER_factors/fastqtl.pc3.peer8.padj.txt \
-output=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2/HCASMC.txt \
-fdr=0.1


# get top eQTLs from GTEx tissue:
bash $scripts/run_get_top_eqtls.gtex.sh \
	346 \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2/


# get top GWAS hits:
Rscript $scripts/get_top_gwas.R \
-input=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt \
-output=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2/GWAS.txt \
-alpha=1e-7


# calculate the LD between GWAS and eQTL hits:
mkdir /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD2
bash $scripts/calculate_ld.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2 \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD2


# plot the number of significant LD variants: 
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2/*.txt > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex2/combined.txt


# run GWAS eQTL LD analysis: 
fdrs=(0.05 0.1 0.2 0.3)
alphas=(1e-8 1e-7 1e-6 1e-5)

for fdr in ${fdrs[@]}; do
for alpha in ${alphas[@]}; do
echo $fdr 
echo $alpha
bash $scripts/run_gwas_eqtl_ld_analysis.sh \
$fdr \
$processed_data/compare_hcasmc_and_gtex_fdr${fdr} \
$alpha \
$processed_data/LD_fdr${fdr}_alpha${alpha}
done 
done 


# run GWAS eQTL LD analysis using D': 
fdrs=(0.05 0.1 0.2 0.3)

for fdr in ${fdrs[@]}; do
echo $fdr 
bash $scripts/run_gwas_eqtl_ld_analysis.gwas_dp0.8.sh \
$fdr \
$processed_data/compare_hcasmc_and_gtex_fdr${fdr} \
/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/CARDIOGRAMplusC4DleadSNPsplusSNPsLD.Dprime0.8 \
$processed_data/LD_fdr${fdr}_dprime0.8
done


# 16/06/18: 
# setup: 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160618
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160618
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160618
mkdir $scripts $processed_data $figures

# format GWAS file for Sherlock:
bash $scripts/format_gwas_for_sherlock.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt \
	$processed_data/sherlock/cad.add.sherlock.txt 


# calculate standard deviation (sd) of gene expression:
Rscript $scripts/calculate_sdY.R \
	-expression=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/combined.filter.norm.rpkm \
	-covariate=../processed_data/160530/find_optimum_num_PEER_factors_matrixeqtl/covariates.matrixeqtl.pc3.peer8.tsv \
	-output=$processed_data/rpkm_sd.pc3.peer8.txt

# run coloc: 
# TBA

# 16/06/24:
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160624
mkdir $scripts


# 16/06/27:
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160627
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160627
mkdir $scripts $processed_data $figures


# run leafcutter:
bash $scripts/run_bam2junc.sh 

python /srv/persistent/bliu2/tools/leafcutter/clustering/leafcutter_cluster.py \
	-j /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments/juncfiles.txt \
	-r /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/leafcutter/ \
	-o leafcutter


# normalize: 
Rscript $scripts/normalize.R \
	-leafcutter=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/leafcutter/leafcutter_perind.counts.gz
	-figure=$figures/intron_missingness_rate.pdf
	-output=$processed_data/leafcutter.norm.tsv


# transpose splicing levelings: 
Rscript /srv/persistent/bliu2/HCASMC_eQTL/scripts/160527/transpose_rpkm.R \
	$processed_data/leafcutter.norm.tsv \
	$processed_data/leafcutter.norm.t.tsv


# find peer factors: 
bash /srv/persistent/bliu2/HCASMC_eQTL/scripts/160527/GetPeerExtended.sh \
	$processed_data/leafcutter.norm.t.tsv \
	$processed_data/


# prepare covariates: 
Rscript /srv/persistent/bliu2/HCASMC_eQTL/scripts/160530/combine_covariates.R \
	--genotype_pc=../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
	--peer=$processed_data/factors.tsv \
	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
	--output=$processed_data/covariates.fastqtl.tsv \
	--gender_coding=letter
bgzip $processed_data/covariates.fastqtl.tsv


# convert splicing level to bed format: 
Rscript $scripts/leafcutter2bed.R \
	-leafcutter=$processed_data/leafcutter.norm.tsv \
	-bed=$processed_data/leafcutter.norm.bed
bgzip $processed_data/leafcutter.norm.bed
tabix -p bed $processed_data/leafcutter.norm.bed.gz


# 16/06/28:
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160628/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160628/
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160628/
mkdir $scripts $processed_data $figures


# filter for biallelic variants:
wd=/srv/persistent/bliu2/HCASMC_eQTL/data/joint3
ln /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_variants.pass.vcf.gz /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_variants.pass.vcf.gz
bcftools view -m2 -M2 -Ov -o $wd/recalibrated_biallelic_variants.pass.vcf $wd/recalibrated_variants.pass.vcf.gz


# beagle imputation without reference: 
bash $scripts/run_beagle_impute.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.pass.vcf \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle


# post-imputation beagle QC:
cat /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.vcf | \
	grep -v "^#" | \
	awk 'BEGIN {FS="\t|;|="; OFS="\t"; print "CHROM","POS","AR2","DR2","AF"} {print $1,$2,$9,$11,$13}' > \
	$processed_data/recalibrated_biallelic_SNP.r2.tsv
Rscript beagle_QC.R \
	-input=$processed_data/recalibrated_biallelic_SNP.r2.tsv \
	-figure_dir=$figures


# update sample names:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
ln /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/old_to_new_sample_name.txt /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/old_to_new_sample_name.txt 
bcftools reheader -s old_to_new_sample_name.txt -o recalibrated_biallelic_variants.beagle.rename.vcf recalibrated_biallelic_variants.beagle.vcf


# filter for variants with dosage R2 >= 0.8:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools view -e 'INFO/DR2<0.8' -o recalibrated_biallelic_variants.beagle.rename.dr2.vcf recalibrated_biallelic_variants.beagle.rename.vcf


# hwe filtering:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
ln /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/caucasian_for_hwe.txt /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/caucasian_for_hwe.txt
grep -v "IMP" recalibrated_biallelic_variants.beagle.rename.dr2.vcf > recalibrated_biallelic_variants.beagle.rename.dr2.2.vcf
vcftools --vcf recalibrated_biallelic_variants.beagle.rename.dr2.2.vcf --keep caucasian_for_hwe.txt --hardy --out hwe_pval
rm recalibrated_biallelic_SNP.beagle.rename.dr2.2.vcf

# select sites with hwe > 1e-6:
tail -n +2 hwe_pval.hwe | awk 'BEGIN{OFS="\t"} {if ($6 >= 1e-6) print $1,$2}' > pass_hwe.txt
bcftools view -T pass_hwe.txt -Ov -o recalibrated_biallelic_variants.beagle.rename.dr2.hwe.vcf recalibrated_biallelic_variants.beagle.rename.dr2.vcf


# filter for MAF > 0.05:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools view \
	--min-af 0.05 --max-af 0.95 \
	-Ov -o recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf \
	recalibrated_biallelic_variants.beagle.rename.dr2.hwe.vcf


# change ID: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools annotate -Oz \
	-o recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	--set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf
tabix -p vcf recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz


# 16/06/29:
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160629/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160629/
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160629/
mkdir $scripts $processed_data $figures


# nominal pass with normalized splicing levels 100kb window: 
bash $scripts/run_fastqtl.nominal.wrap.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.bed.gz \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/covariates.fastqtl.tsv.gz \
	$processed_data/sqtl.nominal.allpairs.normal.1e5.txt \
	"--normal --window 1e5"


# make qqplot and histogram of nominal p-values:
cp 160530/qqplot_pvalue.R $scripts/qqplot_pvalue.R
Rscript $scripts/qqplot_pvalue.R \
	$processed_data/sqtl.nominal.allpairs.txt \
	$figures/qqplot.normal.1e5.pdf \
	$figures/histogram.normal.1e5.pdf


# adjust p-values:
cp 160530/fastqtl_pvalue_corrections.R $scripts/fastqtl_nominal_pvalue_corrections.R
Rscript $scripts/fastqtl_nominal_pvalue_corrections.R $processed_data/sqtl.nominal.allpairs.normal.1e5.txt $processed_data/sqtl.nominal.allpairs.normal.1e5.padj.txt 


# diagnose p-value inflation: 
Rscript /srv/persistent/bliu2/HCASMC_eQTL/scripts/160629/diagnostic_p_value_dist.R


# plot number of significant sqtl vs distance: 
Rscript $scripts/plot_sqtl_vs_distance.R \
	-sqtl_file=$processed_data/sqtl.nominal.allpairs.normal.1e5.padj.txt \
	-all_sqtl_fig=$figures/num_sig_sqtl_vs_dist.pdf \
	-intronic_sqtl_fig=$figures/num_sig_sqtl_within_intron_vs_dist.pdf


# create bedgraph and bigwig files:
parallel -j5 "bash $scripts/160629/bam2bw.sh {}" :::: /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments/sample_list.txt


# merge all bigwig files: 
bash merge_bigwig.sh 


# 16/07/05:
# permutation to calibrate model. 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/
mkdir $scripts/160705 $processed_data/160705 $figures/160705


# permute individual labels in intron level file:
Rscript permute_phenotype.R \
	-input=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160627/leafcutter.norm.bed.gz
	-output=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160705/leafcutter.norm.perm.bed


# run fastqtl on permuted data:
bgzip /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160705/leafcutter.norm.perm.bed
tabix -p bed /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160705/leafcutter.norm.perm.bed.gz
bash $scripts/160629/run_fastqtl.nominal.wrap.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	$processed_data/160705/leafcutter.norm.perm.bed.gz \
	$processed_data/160627/covariates.fastqtl.tsv.gz \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.txt \
	"--normal --window 1e5"

# make qqplot and histogram: 
Rscript $scripts/160629/qqplot_pvalue.R \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.txt \
	$figures/160705/sqtl.perm.qqplot.pdf \
	$figures/160705/sqtl.perm.histogram.pdf
	


# run model with only genotype PCs and genders as covariates:
zcat $processed_data/160627/covariates.fastqtl.tsv.gz | grep -v "InferredCov" | bgzip > $processed_data/160705/covariates.pc.gender.fastqtl.tsv.gz

bash $scripts/160629/run_fastqtl.nominal.wrap.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	$processed_data/160705/leafcutter.norm.perm.bed.gz \
	$processed_data/160705/covariates.pc.gender.fastqtl.tsv.gz \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.noPEER.txt \
	"--normal --window 1e5"

Rscript $scripts/160629/qqplot_pvalue.R \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.noPEER.txt \
	$figures/160705/sqtl.perm.qqplot.noPEER.pdf \
	$figures/160705/sqtl.perm.histogram.noPEER.pdf
	

# run model with PEER factor obtained from top 10000 introns: 
# subset top introns: 
Rscript get_top_introns.R \
	-ratio_file=$processed_data/160627/leafcutter.norm.tsv
	-count_file=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/leafcutter/leafcutter_perind_numers.counts.gz
	-output_file=$processed_data/160705/leafcutter.norm.top10000.tsv

# transpose splicing levelings: 
Rscript /srv/persistent/bliu2/HCASMC_eQTL/scripts/160527/transpose_rpkm.R \
	$processed_data/160705/leafcutter.norm.top10000.tsv \
	$processed_data/160705/leafcutter.norm.top10000.t.tsv

# find peer factors: 
bash /srv/persistent/bliu2/HCASMC_eQTL/scripts/160527/GetPeerExtended.sh \
	$processed_data/160705/leafcutter.norm.top10000.t.tsv \
	$processed_data/160705/

# prepare covariates: 
Rscript /srv/persistent/bliu2/HCASMC_eQTL/scripts/160530/combine_covariates.R \
	--genotype_pc=../processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
	--peer=$processed_data/160705/factors.tsv \
	--sample_info=/srv/persistent/bliu2/HCASMC_eQTL/data/sample_info/sample_info.xlsx \
	--output=$processed_data/160705/covariates.top10000.fastqtl.tsv \
	--gender_coding=letter
bgzip $processed_data/160705/covariates.top10000.fastqtl.tsv

# map sQTLs with PEER factors estimated with top 10000 introns: 
bash $scripts/160629/run_fastqtl.nominal.wrap.sh \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	$processed_data/160705/leafcutter.norm.perm.bed.gz \
	$processed_data/160705/covariates.top10000.fastqtl.tsv.gz \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.top10000.txt \
	"--normal --window 1e5"

# make qqplot and histogram:
Rscript $scripts/160629/qqplot_pvalue.R \
	$processed_data/160705/sqtl.nominal.allpairs.normal.1e5.perm.top10000.txt \
	$figures/160705/sqtl.perm.qqplot.top10000.pdf \
	$figures/160705/sqtl.perm.histogram.top10000.pdf


# map sQTLs with 3,6,9,12,15 PEER factors (obtained from top 10000 introns):
cp $scripts/160530/find_optimal_num_PEER_factors.sh $scripts/160705/find_optimal_num_PEER_factors.sh
bash $scripts/160705/find_optimal_num_PEER_factors.sh


# 160708: 
# run WASP correction
mkdir $scripts/160708 $processed_data/160708 $figures/160708

# copy WASP scripts: 
cp /srv/persistent/bliu2/dase_method/scripts/*WASP* $scripts/160708/


# split multiallelic variants into biallelic records: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/joint3
bcftools norm \
	-f /srv/persistent/bliu2/shared/genomes/hg19/hg19.fa \
	-m -any \
	-Oz -o recalibrated_variants.pass.norm.split.vcf.gz \
	recalibrated_variants.pass.vcf.gz


# generate WASP SNP files: 
python $scripts/160708/generate_WASP_SNP_files.py \
	/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_variants.pass.norm.split.vcf.gz \
	$processed_data/160708/WASP_SNPs/


# find intersecting SNPs: 
python $scripts/160708/WASP_find_intersecting_snps.py \
	/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments/ \
	sample_list.txt


# WASP remap:
# construct table with two columns:
# column 1: genome for pass2 on valk
# column 2: sample directory on durga
bash $scripts/160708/WASP_remap.sh


# WASP filter remapped reads: 
python $scripts/160708/WASP_filter_remapped_reads.py $data/rnaseq2/alignments sample_list.txt 


# clean up: 
bash $scripts/160708/cleanup.sh


# WASP remove duplicate: 
python $scripts/160708/WASP_rmdup.py $data/rnaseq2/wasp $data/rnaseq2/alignments/sample_list.txt 



#---------------- HCASMC specific gene ---------------#
# obj: DE between HCASMC and GTEx
# setup:
mkdir $scripts/160715 $processed_data/160715 $figures/160715


# copy tissue names to HCASMC data directory: 
mkdir $data/gtex
cp /users/joed3/GTExCisEqtls/data/gtex.v6p.eqtl.tissues.txt $data/gtex


# link gtex v6p eQTL data to HCASMC data directory: 
ln -s /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation $data/gtex/v6p
ln -s /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/ $data/gtex


# subsample GTEx tissues for DE:
cp /srv/persistent/bliu2/gtexciseqtls/subsampling/subsample.lists.by.tissue.R $scripts/160715/
mkdir $data/gtex/v6p/subsampling
Rscript $scripts/160715/subsample.lists.by.tissue.bl.R -size=10 # create lists of samples to inlude in the subsets.


# perform subsampling on read counts: 
Rscript $scripts/160715/subsample.R \
	-input='/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_All_Tissues_read_counts_FOR_QC_ONLY.gct' \
	-pattern='*.10.txt' \
	-output_suffix='count' # actually do the subsampling


# perform subsampling on rpkm:
Rscript $scripts/160715/subsample.R \
	-input='/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/v6p/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct' \
	-pattern='*.10.txt' \
	-output_suffix='rpkm'


# merge HCASMC samples into one dataframe: 
bash $scripts/160715/combine_read_counts.sh
bash $scripts/160715/combine_rpkm.sh


# DESeq with one-way ANOVA model:
cp /srv/persistent/bliu2/HCASMC_eQTL/scripts/160614/DSEeq2.fibroblast.R $scripts/160715/find_housekeeping_genes.R 
Rscript $scripts/160715/find_housekeeping_genes.R


# find housekeeping genes by thresholding on mean and variance:
Rscript $scripts/160715/find_housekeeping_genes.2.R


# correct out unwanted variation:
Rscript $scripts/160715/ruvseq.R


# find HCASMC specific genes using the "rank test": 
Rscript $scripts/160715/find_hcasmc_specific_genes.R


# correct out unwanted variation, using size factor corrected counts as input: 
Rscript $scripts/160715/ruvseq.2.R


# find HCASMC specific genes using the "rank test" (without RUVSeq correction):
Rscript $scripts/160715/find_hcasmc_specific_genes.2.R


# HCASMC specific gene biclustering: 
Rscript $scripts/160715/find_hcasmc_specific_genes.3.R



# find HCASMC specific genes using DESeq contrast: 
Rscript $scripts/160715/find_hcasmc_specific_genes.4.R


# GSEA:
cp -r /Users/boshliu/gsea_home/output/aug09/hcasmc_vs_gtex_all.h.all.GseaPreranked.1470779960906/ /Volumes/HCASMC_eQTL/processed_data/160715/hcasmc_vs_gtex_all.h.all
cp -r /Users/boshliu/gsea_home/output/aug09/hcasmc_vs_gtex_all.c2.cp.kegg.GseaPreranked.1470781703555/ /Volumes/HCASMC_eQTL/processed_data/160715/hcasmc_vs_gtex_all.c2.cp.kegg


# Find HCASMC specific gene using information theory: 
bash 160715/prepare_gct_file.sh 
Rscript 160715/find_hcasmc_specific_genes.info_theory.R
Rscript 160715/find_hcasmc_specific_genes.info_theory.quant_norm.R
Rscript 160715/parse_gtf.R
Rscript 160715/manhattan.R

# Overlap HCASMC specific gene and GWAS: 
Rscript 160715/hcasmc_specific_gene_and_GWAS.R 
Rscript 160715/hcasmc_specific_gene_and_GWAS.quant_norm.R 
Rscript 160715/hcasmc_specific_gene_and_GWAS.vary_closest_genes.R 


#### 160724:
# obj: re-map sQTL using WASP files:
# setup: 
mkdir $scripts/160724 $figures/160724 $processed_data/160724

# run leafcutter: 
cp $scripts/160627/run_bam2junc.sh $scripts/160724/

bash $scripts/160724/run_bam2junc.sh

#### 160729:
# run t-SNE visualization
# setup: 
mkdir $scripts/160729 $figures/160729 $processed_data/160729

# create a combined RPKM file of HCASMC and GTEx samples: 
cp $scripts/160603/combine_and_filter_rpkm.R $scripts/160729/combine_and_filter_rpkm.R
Rscript $scripts/160729/combine_and_filter_rpkm.R \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/sample_list.txt \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160729/combined \
	0.1 \
	0.8

# perform t-SNE: 
Rscript $scripts/160729/tsne.2.R


#### 160801: 
# differential expression between serum-fed and serum-starved HCASMC cell lines: 
# setup: 
mkdir $scripts/160801 $figures/160801 $processed_data/160801


# combine read counts: 
cp $scripts/160715/combine_read_counts.sh $scripts/160801/
bash $scripts/160801/combine_read_counts.sh
mv $processed_data/160801/combined.count $processed_data/160801/rnaseq_dase.combined.count
cp /srv/persistent/bliu2/dase/scripts/de_and_db/deseq2.HCASMC_rnaseq_sfFbs_allReps.r \
	$scripts/160801/DESeq2_HCASMC_SF_vs_FBS.R
Rscript $scripts/160801/DESeq2_HCASMC_SF_vs_FBS.R


# get residuals by regression out "batch":
Rscript $scripts/160801/DESeq2_HCASMC_residuals.R


# GSEA enrichment analysis: 
# no script.


#### 160805:
# map eQTLs with fastQTL:
# setup: 
mkdir $scripts/160805 $processed_data/160805 $figures/160805


# find the best combination of genotype PCs and PEER factors:
cp $scripts/160629/run_fastqtl.nominal.wrap.sh $scripts/160805/run_fastqtl.wrap.sh
cp $scripts/160530/find_optimal_num_PEER_factors.sh $scripts/160805
bash $scripts/160805/find_optimal_num_PEER_factors.sh 
cp $scripts/160530/plot_num_egene_vs_cov.R $scripts/160805/
Rscript $scripts/160805/plot_num_egene_vs_cov.R


# create covariate file: 
Rscript $scripts/160530/combine_covariates.R \
	--genotype_pc=$processed_data/160519_genotype_PCA/genotype_pcs.52samples.tsv \
	--peer=$processed_data/160527/factors.tsv \
	--sample_info=$data/sample_info/sample_info.xlsx \
	--output=$processed_data/160805/covariates.pc4.peer8.gender_letter.tsv \
	--gender_coding=letter \
	--num_geno_pc=4 \
	--num_peer_factor=8
bgzip $processed_data/160805/covariates.pc4.peer8.gender_letter.tsv


# nominal pass to map eQTLs: 
bash $scripts/160805/run_fastqtl.wrap.sh \
	$data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	$processed_data/160530/combined.filter.norm.bed.gz \
	$processed_data/160805/covariates.pc4.peer8.gender_letter.tsv.gz \
	$processed_data/160805/hcasmc.eqtl.pc4.peer8.txt \
	"--window 1e6"

# select variants with p-value less than 1e-3: 
zcat $processed_data/160805/hcasmc.eqtl.pc4.peer8.txt.gz | awk 'BEGIN{OFS="\t";print "snp","pval"}{if ($4<1e-3) {print $2,$4}}' >  $processed_data/160805/hcasmc.eqtl.pc4.peer8.pval1e-3.txt


# run p-value correction:
Rscript $scripts/160629/fastqtl_nominal_pvalue_corrections.R $processed_data/160805/hcasmc.eqtl.pc4.peer8.txt $processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt


# plot the p-values: 
Rscript $scripts/160629/qqplot_pvalue.R \
	$processed_data/160805/hcasmc.eqtl.pc4.peer8.txt \
	$figures/160805/hcasmc.eqtl.qqplot.pdf \
	$figures/160805/hcasmc.eqtl.histogram.pdf


# plot the number of eQTLs vs distance to TSS: 
cp $scripts/160629/plot_sqtl_vs_distance.R $scripts/160805/plot_eqtl_vs_distance.R
Rscript $scripts/160805/plot_eqtl_vs_distance.R \
	-eqtl_file=$processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt \
	-num_eqtl_vs_dist=$figures/160805/num_sig_eqtl_vs_dist.pdf \
	-pval_vs_dist=$figures/160805/eqtl_pval_vs_dist.pdf


# permutation pass to map eQTLs: 
bash $scripts/160805/run_fastqtl.wrap.sh \
	$data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf.id.gz \
	$processed_data/160530/combined.filter.norm.bed.gz \
	$processed_data/160805/covariates.pc4.peer8.gender_letter.tsv.gz \
	$processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.txt \
	"--window 1e6 --permute 1000 10000"


# correct permutation eQTL pvalue:
cp $scripts/160530/fastqtl_pvalue_corrections.R $scripts/160805/fastqtl_permuted_pvalue_corrections.R
Rscript $scripts/160805/fastqtl_permuted_pvalue_corrections.R \
	../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.txt \
	../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.padj.txt


# create symbolic link of eQTL data: 
ln -s /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_allpairs_FOR_QC_ONLY/ /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/
cp /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/Metasoft_tissue_order.txt /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160805/
# manually added HCASMC as the first line in Metasoft_tissue_order.txt


# reformat the genotype column of HCASMC eqtl file: 
gzip $processed_data/160805/hcasmc.eqtl.pc4.peer8.txt
Rscript $scripts/160805/format_hcascm_eqtl.R # output $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt


#----------------- Metasoft -----------------
## Run Metasoft (full sample, select eQTL with p<1e-3): 
# Copy HCASMC eQTL data to appropriate location: 
gzip $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt
ln -s $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt.gz $processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/HCASMC_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz

# Concatenate eQTL result for all tissue, append tissue name, and sort according to gene ID and SNP: 
bash $scripts/160805/concatenate_eqtl_tissues.sh 

# Make input file for Metasoft:
mkdir $processed_data/160805/metasoft_input/
cat $processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py > $processed_data/160805/metasoft_input/metasoft_input.txt

# Split metasoft input by chromosome (to reduce memory footprint for next step): 
bash $scripts/160805/split_metasoft_input_by_chr.sh

# Run METASOFT:
mkdir $processed_data/160805/metasoft_output/
parallel -j12 bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input/metasoft_input.{}.txt $processed_data/160805/metasoft_output/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output/metasoft_output.{}.mcmc.log ::: {1..22} X

# merge metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt


## Run Metasoft (subsampled to 52, select eQTL with p<1e-3): 
# copy HCASMC eQTL data to appropriate location:
mkdir $processed_data/160816/subsampling/HCASMC
ln -s $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt.gz $processed_data/160816/subsampling/HCASMC/HCASMC_52.allpairs.txt.gz

# Concatenate eQTL result for all tissue, append tissue name, and sort according to gene ID and SNP: 
bash $scripts/160805/concatenate_eqtl_tissues.subsample.sh 

# Generate Metasoft input file:
mkdir $processed_data/160805/metasoft_input_subsample_52/ 
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py > $processed_data/160805/metasoft_input_subsample_52/metasoft_input.txt

# Split metasoft input by chromosome: 
parallel 'grep "_{}_" ../processed_data/160805/metasoft_input_subsample_52/metasoft_input.txt > ../processed_data/160805/metasoft_input_subsample_52/metasoft_input.{}.txt' ::: {1..22} X

# Run METASOFT:
mkdir $processed_data/160805/metasoft_output_subsample_52
parallel bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input_subsample_52/metasoft_input.{}.txt $processed_data/160805/metasoft_output_subsample_52/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output_subsample_52/metasoft_output.{}.mcmc.log ::: {1..22} X

# Merge Metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt

## Run Metasoft (subsample to 52, select eQTL with p<1e-2)
# generate metasoft input file:
mkdir $processed_data/160805/metasoft_input_subsample_52_p1e-2/
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py 1e-2 > $processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.txt

# split Metasoft input by chromosome: 
parallel 'grep "_{}_" ../processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.txt > ../processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.{}.txt' ::: {1..22} X

# run METASOFT:
mkdir $processed_data/160805/metasoft_output_subsample_52_p1e-2
parallel bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.{}.txt $processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.{}.mcmc.log ::: {1..22} X

# Merge Metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt

# plot heatmap of m-values: 
Rscript $scripts/160805/plot_mvalue_heatmap.R 


#### 160811:
# fine maping by inspection: 
# setup: 
mkdir $scripts/160811 $processed_data/160811 $figures/160811


# extract known gwas hits from the supplemetary excel of Nikpay 2015 NG: 
Rscript $scripts/160811/extract_gwas_loci.R  # ../processed_data/160811/gwas_loci.cad.all.FWER.txt

# format the GWAS loci file so it contains one column chr_pos:
cat ../processed_data/160811/gwas_loci.cad.all.FWER.txt | awk 'NR>1 {print $1"_"$2}' > ../processed_data/160811/gwas_loci.cad.all.FWER.chr_pos.txt 

# select top GWAS variants from eQTL:
grep -f ../processed_data/160811/gwas_loci.cad.all.FWER.chr_pos.txt ../processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt > ../processed_data/160811/tested_genes_at_gwas_top_hits.txt

# plot p-values of gene:snp pairs at GWAS top hits: 
Rscript $scripts/160811/plot_gwas_hits_eqtl_pval.R

# subset to eGene with p-value < 0.05 at GWAS top hits: 
cat ../processed_data/160811/tested_genes_at_gwas_top_hits.txt | awk '{if ($4<0.05) print $1}' > ../processed_data/160811/egenes_at_gwas_top_hits.p05.txt

# subset to eQTL file to eGenes selected above: 
grep -f ../processed_data/160811/egenes_at_gwas_top_hits.p05.txt ../processed_data/160805/hcasmc.eqtl.pc4.peer8.padj.txt | sort -k1,2 -V > ../processed_data/160811/tested_snps_at_egenes.txt

# add rs ID to each snp: 
cat ../processed_data/160811/tested_snps_at_egenes.txt | awk 'BEGIN {FS="\t|_";OFS="\t"} {print $2,$3}' > ../processed_data/160811/tested_snps_at_egenes.chr_pos.txt
cat ../processed_data/160811/tested_snps_at_egenes.chr_pos.txt | python $scripts/160811/subset_dbsnp.py > ../processed_data/160811/tested_snps_at_egenes.dbsnp146.txt
Rscript $scripts/160811/add_rsid.R -dbsnp=../processed_data/160811/tested_snps_at_egenes.dbsnp146.txt -eqtl=../processed_data/160811/tested_snps_at_egenes.txt -out=../processed_data/160811/tested_snps_at_egenes.rsid.txt


# make locuszoom plot for ENSG00000197208.5 at rs273909 (SLC22A4-SLC22A5):
grep ENSG00000197208.5 ../processed_data/160811/tested_snps_at_egenes.rsid.txt | awk 'BEGIN{OFS="\t"; print "markername","pval"} {print $7,$4}'> ../processed_data/160811/ENSG00000197208.metal.txt
echo -e "snp\tstring\tcolor" > ../processed_data/160811/ENSG00000197208.markers.txt; echo -e "rs1045020\teQTL\tblue" >> ../processed_data/160811/ENSG00000197208.markers.txt; echo -e "rs273909\tGWAS\tpink" >> ../processed_data/160811/ENSG00000197208.markers.txt
locuszoom --metal ../processed_data/160811/ENSG00000197208.metal.txt --pvalcol pval --markercol markername --refsnp rs273909 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR --denote-markers-file ../processed_data/160811/ENSG00000197208.markers.txt title="eQTL (ENSG00000197208,SLC22A4)" --prefix ../processed_data/160811/eqtl_ENSG00000197208_rs273909


# cast CAD GWAS into METAL format: 
cat $data/gwas/cad.add.160614.website.txt | awk 'BEGIN {OFS="\t"} {print $1, $11}' > $data/gwas/cad.add.160614.website.metal.txt
grep -f ../processed_data/160811/tested_snps_at_egenes.chr_pos.txt  /srv/persistent/bliu2/shared/dbsnp/snp146.txt > ../processed_data/160811/tested_snps_at_egenes.dbsnp146.txt


# make locuszoom plot for GWAS hit rs273909 (SLC22A4-SLC22A5):
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs273909 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR --denote-markers-file ../processed_data/160811/ENSG00000197208.markers.txt title="GWAS (rs273909)" --prefix ../processed_data/160811/gwas_rs273909_ENSG00000197208


# make forest PM plot for rs273909:
cd /srv/persistent/bliu2/tools/ForestPMPlot/
python /srv/persistent/bliu2/tools/ForestPMPlot/pmplot.py \
$processed_data/160805/metasoft_input/metasoft_input.5.txt \
$processed_data/160805/metasoft_output/metasoft_output.5.mcmc.txt \
$processed_data//160805/Metasoft_tissue_order.alphabetical.txt \
$processed_data/160805/Metasoft_tissue_idx.txt \
ENSG00000197208.5_5_131667353_A_G_b37 \
ENSG00000197208,SLC22A4 \
$figures/160811/rs273909_ENSG00000197208.pdf

# make locuszoom plot for eQTL ENSG00000182511.7 at rs2521501 (SLC22A4-SLC22A5):
grep ENSG00000182511.7 ../processed_data/160811/tested_snps_at_egenes.rsid.txt | awk 'BEGIN{OFS="\t"; print "markername","pval"} {print $7,$4}'> ../processed_data/160811/ENSG00000182511.metal.txt
echo -e "snp\tstring\tcolor" > ../processed_data/160811/ENSG00000182511.markers.txt; echo -e "rs1045020\teQTL\tblue" >> ../processed_data/160811/ENSG00000182511.markers.txt; echo -e "rs273909\tGWAS\tpink" >> ../processed_data/160811/ENSG00000182511.markers.txt
locuszoom --metal ../processed_data/160811/ENSG00000182511.metal.txt --pvalcol pval --markercol markername --refsnp rs2521501 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000182511,FES)" --prefix ../processed_data/160811/eqtl_ENSG00000182511_rs2521501


# make locuszoom plot for GWAS hit rs2521501 (FURIN-FES):
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs2521501 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="GWAS (rs2521501)" --prefix ../processed_data/160811/gwas_ENSG00000182511_rs2521501


# make forest PM plot for rs2521501: 
cd /srv/persistent/bliu2/tools/ForestPMPlot/
python /srv/persistent/bliu2/tools/ForestPMPlot/pmplot.py \
$processed_data/160805/metasoft_input/metasoft_input.15.txt \
$processed_data/160805/metasoft_output/metasoft_output.15.mcmc.txt \
$processed_data//160805/Metasoft_tissue_order.alphabetical.txt \
$processed_data/160805/Metasoft_tissue_idx.txt \
ENSG00000182511.7_15_91437388_A_T_b37 \
ENSG00000182511,FES \
$figures/160811/rs2521501_ENSG00000182511.pdf


#### 160813
#### convert vcf to bed fasta pairs
# setup:
mkdir $scripts/160813

# generate sample list: 
python $scripts/160813/gen_sample_list.py > $data/joint3/sample_list.txt

# convert vcf to bed fasta pairs: 
parallel -a $data/joint3/sample_list.txt -j20 \
	python $scripts/160813/vcf2bedfa.py \
	$data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.vcf \
	{} \
	$data/joint3/bed_fa/{}


#### 160816
#### eQTL mapping on subsampled GTEx tissues:
# setup: 
mkdir $scripts/160816 $processed_data/160816 $figures/160816


# copy fastqtl subsampling code: 
cp /srv/persistent/bliu2/gtexciseqtls/subsampling/nominal.pass.subsamples.v6p.{core.sh,sh} $scripts/160816
cp /srv/persistent/bliu2/gtexciseqtls/subsampling/run_FastQTL_threaded_subsample.py  $scripts/160816
cp /srv/persistent/bliu2/gtexciseqtls/subsampling/subsample.lists.by.tissue.R $scripts/160816


# subset to covariates to first 10: 
bash $scripts/160816/subset_covariates.sh 


# subsample tissues: 
Rscript $scripts/160816/subsample.lists.by.tissue.R


# run fastqtl nominal pass on subsampled GTEx tissues: 
bash $scripts/160816/nominal.pass.subsamples.v6p.sh


# run fastqtl permutation pass on subsampled GTEx tissues: 
bash $scripts/160816/permutation.pass.subsamples.v6p.sh


#### replication of HCASMC eQTLs in GTEx tissues: 
# select HCASMC eQTLs with FDR < 0.05: 
cat ../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.padj.txt | \
awk 'BEGIN{OFS="\t"} {if ($19<=0.05) {print $1,$6"_b37"}}' | \
sed "s/chr//" > ../processed_data/160816/hcasmc.eqtl.pc4.peer8.perm.padj.fdr05.txt


# subset GTEx tissue to association significant in HCASMC:
mkdir $processed_data/160816/replication/
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
python $scripts/160816/subset_eQTLs.py \
../processed_data/160816/hcasmc.eqtl.pc4.peer8.perm.padj.fdr05.txt \
$data/gtex/v6p/v6p_fastQTL_allpairs_FOR_QC_ONLY/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz \
../processed_data/160816/replication/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.sig_in_HCASMC.txt


# calculate pi1 statistics for GTEx tissues
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
Rscript $scripts/160816/pi1_calc.R \
../processed_data/160816/replication/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.sig_in_HCASMC.txt \
{} > ../processed_data/160816/replication/pi1.txt


# plot pi1 statistic:
Rscript $scripts/160816/pi1_plot.R \
../processed_data/160816/replication/pi1.txt \
../figures/160816/pi1.pdf


# subset subsampled GTEx tissue to association significant in HCASMC:
mkdir $processed_data/160816/replication_subsample/
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
python $scripts/160816/subset_eQTLs.py \
../processed_data/160816/hcasmc.eqtl.pc4.peer8.perm.padj.fdr05.txt \
$processed_data/160816/subsampling/{}/{}_52.allpairs.txt.gz \
../processed_data/160816/replication_subsample/{}_52.allpairs.sig_in_HCASMC.txt


# calculate pi1 statistics for GTEx tissues
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
Rscript $scripts/160816/pi1_calc.R \
../processed_data/160816/replication_subsample/{}_52.allpairs.sig_in_HCASMC.txt \
{} > ../processed_data/160816/replication_subsample/pi1.txt


# plot pi1 statistic:
Rscript $scripts/160816/pi1_plot.R \
../processed_data/160816/replication_subsample/pi1.txt \
../figures/160816/pi1_subsample.pdf


#### 160824
#### run eCavier
# setup: 
mkdir $scripts/160824 $processed_data/160824 $figures/160824

# check eQTL file is sorted: 
zcat ../processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt.gz | sort --parallel=10 -k1,1 -k2,2 -V | gzip > ../processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.sort.txt.gz


# split eqtl file by gene name: 
mkdir ../processed_data/160824/eqtl_by_gene/
zcat ../processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.sort.txt.gz | python $scripts/160824/split_eqtl_by_gene.py ../processed_data/160824/eqtl_by_gene/
for input in $(ls ../processed_data/160824/eqtl_by_gene/ENSG*txt); do
	cat $input | awk '{print $2,$5/$6}' > ${input/txt/eqtl.zscore}
done

# change rsid to chr_pos_ref_alt_b37 in all 1000 Genomes p1v3 vcf files:
bash $scripts/160824/reformat_rsid.sh
parallel -j11 tabix -p vcf /srv/persistent/bliu2/shared/1000genomes/phase1v3/ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.rsid.vcf.gz ::: {1..22}

# rsids for indels are labeled by chr:pos:<D|I>; change these rsids to 
# be consisten with 1000 Genomes phase 1 version 3: 
# cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt | python $scripts/160824/update_gwas_rsid.py > /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.rsid.txt

# intersect gwas and each eqtl (split by gene):
Rscript $scripts/160824/intersect_gwas_eqtl.R


# subset to genes with at least 1 eQTL and 1 GWAS hit (p<1e-3):
cat ../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.txt | awk '{print $1}' > ../processed_data/160824/gene_id.txt
parallel -a ../processed_data/160824/gene_id.txt -j40 Rscript \
	$scripts/160824/calc_min_pval.R \
	../processed_data/160824/eCAVIAR_input/{}.eqtl.zscore \
	../processed_data/160824/eCAVIAR_input/{}.gwas.zscore \
	{} > ../processed_data/160824/gene_min_pval.txt
Rscript $scripts/160824/select_gene_by_eqtl_and_gwas_pval.R > ../processed_data/160824/gene_id_gwas1e-4_eqtl1e-4.txt


# for each <gene_id>.gwas.zscore: 
# 	calculate LD (using plink and European ancestry)

# run caviar: 
# eqtl file: /srv/persistent/bliu2/HCASMC_eQTL/processed_data//160805/hcasmc.eqtl.pc4.peer8.padj.txt
# cad file: /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt


#### 160829
#### eCAVIAR
# setup: 
mkdir $scripts/eCAVIAR $processed_data/eCAVIAR $figures/eCAVIAR


# link GWAS hits: 
ln ../processed_data/160811/gwas_loci.cad.all.FWER.txt $processed_data/eCAVIAR/gwas_loci.cad.all.FWER.txt


# select GWAS loci (variants less than 50 variants away from the GWAS top hits): 
Rscript $scripts/eCAVIAR/select_gwas_loci.R


# select eGenes at GWAS loci:
Rscript $scripts/eCAVIAR/select_eGene.R


# run eCAVIAR:
bash $scripts/eCAVIAR/ecaviar.sh


# make locuszoom plot:
parallel Rscript $scripts/eCAVIAR/plot_locus.R \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data/eCAVIAR/eCAVIAR_input/{}.gwas.zscore \
	/srv/persistent/bliu2/HCASMC_eQTL/processed_data//eCAVIAR/eCAVIAR_input/{}.eqtl.zscore \
	../figures/eCAVIAR/{}.pdf ::: ENSG00000100014.15 ENSG00000197208.5 ENSG00000198270.8 ENSG00000226972.2 ENSG00000234380.1 ENSG00000236838.2


# make locuszoom plot: 
parallel --xapply $scripts/eCAVIAR/locuszoom.sh {1} {2} $figures/eCAVIAR/ ::: ENSG00000100014.15 ENSG00000197208.5 ENSG00000198270.8 ENSG00000226972.2 ENSG00000234380.1 ENSG00000236838.2 ::: rs5760295 rs273909 rs77684561 rs4245791 rs8134775 rs216172


#### eCAVIAR for downsampled GTEx tissue:

# extract top gwas hits from the supplemetary excel of Nikpay 2015 NG: 
cp $scripts/160811/extract_gwas_loci.R $scripts/eCAVIAR/extract_gwas_variants.R
Rscript $scripts/eCAVIAR/extract_gwas_variants.R \
	/srv/persistent/bliu2/HCASMC_eQTL/data/gwas/nikpay_2015_ng.xlsx \
	../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt


# select GWAS loci: 
Rscript $scripts/eCAVIAR/select_gwas_loci.R \
	../processed_data/eCAVIAR/gwas_loci.cad.all.genomewide_fdr_merged.txt \
	../processed_data/eCAVIAR/cad_gwas_loci/ \
	50 


# select GWAS loci (CAD + MI + recessive + metabochip + all known): 
Rscript $scripts/eCAVIAR/select_gwas_loci.R \
	../data/gwas/gwas.cad.mi.recessive.metabochip.txt \
	../data/gwas/cad_gwas_loci_combined_w100/ \
	100

Rscript $scripts/eCAVIAR/select_gwas_loci.R \
	../data/gwas/gwas.cad.mi.recessive.metabochip.txt \
	../data/gwas/cad_gwas_loci_combined_w50/ \
	50


# preprocessing for selecting eGenes: 
parallel -j20 -a $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt \
	bash $scripts/eCAVIAR/select_eGene.preprocess.sh \
	../processed_data/160816/subsampling/{}/{}_52.allpairs.txt.gz \
	../processed_data/160816/subsampling/{}/{}_52.allpairs.sid_parsed.txt


# split gene-snp pairs by chromsomes:
parallel -j11 bash $scripts/eCAVIAR/split_eqtl_by_chr.sh \
	../processed_data/160816/subsampling/{1}/{1}_52.allpairs.sid_parsed.txt \
	../processed_data/160816/subsampling/{1}/{1}_52.allpairs.sid_parsed.{2}.txt \
	{2} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt ::: {1..22}


# select eGenes:
mkdir ../processed_data/eCAVIAR/eCAVIAR_input4
parallel -j5 Rscript $scripts/eCAVIAR/select_eGene.separate_gwas_loci.R \
	../processed_data/160816/subsampling/{1}/{1}_52.allpairs.sid_parsed.{2}.txt \
	../processed_data/eCAVIAR/eCAVIAR_input4/{1}/ \
	{2} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt ::: {1..22}


# run eCAVIAR:
mkdir ../processed_data/eCAVIAR/eCAVIAR_output4
parallel -j15 bash $scripts/eCAVIAR/ecaviar.sh \
	../processed_data/eCAVIAR/eCAVIAR_input4/{} \
	../processed_data/eCAVIAR/eCAVIAR_output4/{} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt


# count eCAVIAR hits: 
subl $scripts/eCAVIAR/analyze_ecaviar_result.sh 


# make locuszoom plot: 
bash $scripts/eCAVIAR/locuszoom.sh ENSG00000118526.6 rs2327429 $figures/eCAVIAR/
parallel --xapply $scripts/eCAVIAR/locuszoom.sh {1} {2} $figures/eCAVIAR/ ::: ENSG00000118526.6 ENSG00000118526.6 ENSG00000188735.8 ENSG00000198270.8 ENSG00000198270.8 ENSG00000226972.2 ENSG00000234380.1 ENSG00000236838.2 ENSG00000257218.1 ::: rs6569913 rs2327429 rs148608463 rs76741465 rs77684561 rs539702042 rs8134775 rs216172 rs2464190

Rscript $scripts/eCAVIAR/locuszoom.preprocessing.R
parallel --xapply $scripts/eCAVIAR/locuszoom.sh {1} {2} $figures/eCAVIAR/ :::: <(cut -f4,4 ../processed_data/eCAVIAR/locuszoom/HCASMC.txt) :::: <(cut -f3,3 ../processed_data/eCAVIAR/locuszoom/HCASMC.txt)


#### end eCAVIAR for downsampled GTEx tissue: 


#### eCAVIAR for full-sample GTEx tissue:

# preprocessing for selecting eGenes: 
parallel -j11 -a $data/gtex/gtex.v6p.eqtl.tissues.txt \
	bash $scripts/eCAVIAR/select_eGene.preprocess.sh \
	../data/gtex/v6p/v6p_fastQTL_allpairs_FOR_QC_ONLY/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz \
	../processed_data/160816/fullsample/{}/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.sid_parsed.txt


# split gene-snp pairs by chromsomes:
parallel -j11 bash $scripts/eCAVIAR/split_eqtl_by_chr.sh \
	../processed_data/160816/fullsample/{1}/{1}_Analysis.v6p.FOR_QC_ONLY.allpairs.sid_parsed.txt \
	../processed_data/160816/fullsample/{1}/{1}_Analysis.v6p.FOR_QC_ONLY.allpairs.sid_parsed.{2}.txt \
	{2} :::: $data/gtex/gtex.v6p.eqtl.tissues.txt ::: {1..22}


# select eGenes:
# mkdir ../processed_data/eCAVIAR/eCAVIAR_input_fullsample
parallel -j11 Rscript $scripts/eCAVIAR/select_eGene.separate_gwas_loci.R \
	../processed_data/160816/fullsample/{1}/{1}_Analysis.v6p.FOR_QC_ONLY.allpairs.sid_parsed.{2}.txt \
	../processed_data/eCAVIAR/eCAVIAR_input_fullsample/{1}/ \
	{2} :::: $data/gtex/gtex.v6p.eqtl.tissues.txt ::: {1..22}


# run eCAVIAR:
# mkdir ../processed_data/eCAVIAR/eCAVIAR_output_fullsample
parallel -j11 bash $scripts/eCAVIAR/ecaviar.sh \
	../processed_data/eCAVIAR/eCAVIAR_input_fullsample/{} \
	../processed_data/eCAVIAR/eCAVIAR_output_fullsample/{} :::: $data/gtex/gtex.v6p.eqtl.tissues.txt


# # count eCAVIAR hits: 
# subl $scripts/eCAVIAR/analyze_ecaviar_result.sh 


# # make locuszoom plot: 
# bash $scripts/eCAVIAR/locuszoom.sh ENSG00000118526.6 rs2327429 $figures/eCAVIAR/
# parallel --xapply $scripts/eCAVIAR/locuszoom.sh {1} {2} $figures/eCAVIAR/ ::: ENSG00000118526.6 ENSG00000118526.6 ENSG00000188735.8 ENSG00000198270.8 ENSG00000198270.8 ENSG00000226972.2 ENSG00000234380.1 ENSG00000236838.2 ENSG00000257218.1 ::: rs6569913 rs2327429 rs148608463 rs76741465 rs77684561 rs539702042 rs8134775 rs216172 rs2464190

# Rscript $scripts/eCAVIAR/locuszoom.preprocessing.R
# parallel --xapply $scripts/eCAVIAR/locuszoom.sh {1} {2} $figures/eCAVIAR/ :::: <(cut -f4,4 ../processed_data/eCAVIAR/locuszoom/HCASMC.txt) :::: <(cut -f3,3 ../processed_data/eCAVIAR/locuszoom/HCASMC.txt)

#### end eCAVIAR for full-sample GTEx tissue

#### eCAVIAR for HCASMC eQTL mapped with RASQUAL: 
# split eQTL by chromosome: 
parallel -j11 bash $scripts/eCAVIAR/split_eqtl_by_chr.sh \
	../processed_data/rasqual/merged_output/prioritized_genes.allpairs.txt \
	../processed_data/rasqual/merged_output/prioritized_genes.allpairs.{}.txt \
	{} ::: {1..22}

# use 1000G LD for both GWAS and eQTL:
parallel Rscript $scripts/eCAVIAR/select_eGene.separate_gwas_loci.rasqual.R \
	../processed_data/rasqual/merged_output/prioritized_genes.allpairs.{}.txt \
	../data/gwas/cad_gwas_loci_combined_w50 \
	../processed_data/eCAVIAR/eCAVIAR_input_rasqual_110716/ \
	{} ::: {1..22}

bash $scripts/eCAVIAR/ecaviar_bak.sh \
	../processed_data/eCAVIAR/eCAVIAR_input_rasqual_110716/ \
	../processed_data/eCAVIAR/eCAVIAR_output_rasqual_110716/

# use 1000G LD for only GWAS (100 SNP window): 
parallel -j10 Rscript $scripts/eCAVIAR/select_eGene.separate_gwas_loci.rasqual.R \
	../processed_data/rasqual/merged_output/prioritized_genes.allpairs.{}.txt \
	../data/gwas/cad_gwas_loci_combined_w100 \
	../processed_data/eCAVIAR/eCAVIAR_input_rasqual_111516/ \
	{} ::: {1..22}

bash $scripts/eCAVIAR/ecaviar.sh \
	../processed_data/eCAVIAR/eCAVIAR_input_rasqual_111516/ \
	../processed_data/eCAVIAR/eCAVIAR_output_rasqual_111516/ 


# use 1000G LD for only GWAS (50 SNP window): 
parallel Rscript $scripts/eCAVIAR/select_eGene.separate_gwas_loci.rasqual.R \
	../processed_data/rasqual/merged_output/prioritized_genes.allpairs.{}.txt \
	../data/gwas/cad_gwas_loci_combined_w50 \
	../processed_data/eCAVIAR/eCAVIAR_input_rasqual_111616/ \
	{} ::: {1..22}

bash $scripts/eCAVIAR/ecaviar.sh \
	../processed_data/eCAVIAR/eCAVIAR_input_rasqual_111616/ \
	../processed_data/eCAVIAR/eCAVIAR_output_rasqual_111616/ 


# ls *col | sed 's/\.ecaviar_col//' > sample_list.txt
parallel -j1 bash eCAVIAR/ecaviar2locuszoom.sh \
	../processed_data/eCAVIAR/eCAVIAR_output_rasqual/{}.ecaviar_col \
	../data/gwas/cad.add.160614.website.txt \
	../processed_data/rasqual/output/{}.pval.txt \
	../processed_data/eCAVIAR/tmp \
	../figures/locuszoom :::: ../processed_data/eCAVIAR/eCAVIAR_output_rasqual/sample_list.txt




#### eGenes vs sample size:

# setup: 
mkdir $scripts/egenes_vs_sample_size $figures/egenes_vs_sample_size $processed_data/egenes_vs_sample_size

# copy some scripts: 
cp $scripts/160530/plot_num_egenes_vs_sample_size.R $scripts/egenes_vs_sample_size


# plot the number of egenes vs sample size: 
Rscript $scripts/egenes_vs_sample_size/plot_num_egenes_vs_sample_size.R


# count the number of egenes in subsampled (52 individual) GTEx tissues:
bash $scripts/egenes_vs_sample_size/count_num_egenes_subsampled_52.sh


# plot the number of egenes in subsampled tissues:
Rscript $scripts/egenes_vs_sample_size/plot_num_egenes_subsampled_52.R

#### end eGenes vs sample size

#---------- HCASMC specific eQTLs -------------#
# setup:
mkdir hcasmc_specific_eqtl ../processed_data/hcasmc_specific_eqtl ../figures/hcasmc_specific_eqtl


# find HCASMC specific eQTLs (full sample): 
parallel -j5 Rscript ../scripts/hcasmc_specific_eqtl/find_tissue_specific_eqtls.R ../processed_data/160805/metasoft_output/metasoft_output.{2}.mcmc.txt ../figures/hcasmc_specific_eqtl/fullsample/{1}/ ../processed_data/hcasmc_specific_eqtl/fullsample/{1}/stat.{2}.txt ../processed_data/hcasmc_specific_eqtl/fullsample/{1}/tissue_specific_eqtl.{2}.txt {1} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt ::: {1..22}


# find tisue specific eQTLs (subsampled to 52):
parallel -j5 Rscript ../scripts/hcasmc_specific_eqtl/find_tissue_specific_eqtls.R ../processed_data/160805/metasoft_input_subsample_52/metasoft_output.{2}.mcmc.txt ../figures/hcasmc_specific_eqtl/subsample_52/{1}/ ../processed_data/hcasmc_specific_eqtl/subsample_52/{1}/stat.{2}.txt ../processed_data/hcasmc_specific_eqtl/subsample_52/{1}/tissue_specific_eqtl.{2}.txt {1} :::: $data/gtex/gtex.v6p.eqtl.tissues.with_hcasmc.txt ::: {1..22}


# Plot examples of HCASMC-specific eQTL:
bash hcasmc_specific_eqtl/plot_examples_of_hcasmc_eqtl.sh


# make manhattan plot: 
Rscript hcasmc_specific_eqtl/manhattan.R


# Calculate eQTL specificity using information theory (fill NAs with 1e-8):
cat ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python $scripts/hcasmc_specific_eqtl/calc_eqtl_specificity.py ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/eqtl.1e-8.txt.gz ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/specificity.1e-8.txt.gz 1e-8


# Calculate eQTL specificity using information theory (fill NAs using mean imputation):
cat ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python $scripts/hcasmc_specific_eqtl/calc_eqtl_specificity.py ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/eqtl.mean.txt.gz ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/specificity.mean.txt.gz -1


# select HCASMC specific eQTL based on QSI, and plot distribution of their p-values: 
bash hcasmc_specific_eqtl/plot_hcasmc_specific_eqtls_pval.sh


# naively select HCASMC specific eQTL: 
cat ../processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python hcasmc_specific_eqtl/naive_select.py > ../processed_data/hcasmc_specific_eqtl/hcasmc_specific_eqtl.naive_select.txt


# Make forest plot: 
bash hcasmc_specific_eqtl/plot_hcasmc_specific_eqtls.sh 




#------------- ATACseq --------------#
# setup: 
mkdir atacseq ../processed_data/atacseq ../figures/atacseq


# prepare data:
bash atacseq/transfer.sh


# process all FBS ATACseq data with Kundaje pipeline (need to run on nandi):
bash atacseq/run_fbs.sh
bash atacseq/run_sf.sh


# Run WASP:
cd atacseq/wasp/; snakemake --cores 10


# Count ASE (ASVCF):
bash atacseq/asvcf/vcf2asvcf.sh
bash atacseq/asvcf/concordance.sh


# Count peak coverage: 
bash atacseq/count/define_peak.sh


# RASQUAL:
bash atacseq/rasqual/make_input.sh
bash atacseq/rasqual/rasqual.wrapper.sh



#------------ GWAS ATACseq overlap -------------#
# setup: 
mkdir gwas_atacseq_overlap ../processed_data/gwas_atacseq_overlap ../figures/gwas_atacseq_overlap


# Method 1: Overlap with all GWAS variants passing threshold:
# plot fraction of overlaps: 
Rscript gwas_atacseq_overlap/gwas_thresholding/gwas_thresholding.R


# Method 2: GREGOR
# Preprocessing: 
Rscript gwas_atacseq_overlap/gregor/LD_proxy.R 
Rscript gwas_atacseq_overlap/gregor/merge_peaks.R
python gwas_atacseq_overlap/gregor/parse_encode_vocabulary.py > ../data/encode/dnase_seq_2007_2012/controlled_vocabulary/human_cell_types.txt

# Perform overlap and calculate enrichment statistics: 
Rscript gwas_atacseq_overlap/gregor/overlap_enrichment.R


#------------ GWAS eQTL overlap ---------------#
mkdir gwas_eqtl_overlap ../processed_data/gwas_eqtl_overlap ../figures/gwas_eqtl_overlap


# Method 1: Overlap with all GWAS variants pass threshold.
Rscript gwas_eqtl_overlap/gwas_thresholding/overlap.metasoft.R


# Method 2: Gregor:
Rscript gregor/overlap_enrichment.R

#### end GWAS eQTL overlap


#### average RNAseq read depth: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments
grep "Number of input reads" */Log.final.out | awk 'BEGIN{FS="|"}{print $2}' | awk '{print $1}' > read_depths.txt 
#### end average RNAseq read depth


#### average ATACseq read depth: 
cd /users/bliu2/atacseq/2305/out/qc
grep "reads" */*align.log | awk 'BEGIN{FS=":| reads"}{print $2}' 
#### end average ATACseq read depth


#### tarid: 
mkdir tarid ../processed_data/tarid ../figures/tarid


# get rpkm for TARID and TCF21:
bash tarid/tcf21_tarid_expression.sh 


# get sample ID to tissue table:
bash tarid/gen_sampleID_to_tissue_table.sh


# plot TARID and TCF21 expresssion:
Rscript tarid/tarid_vs_tcf21_expression.R # ../figures/tarid_vs_tcf21.pdf, ../figures/tarid_vs_tcf21.median.pdf


# make PM plot for TCF21 and rs2327429 (full sample): 
grep "ENSG00000118526.6_6_134209837_T_C_b37" ../processed_data/160805/metasoft_output/metasoft_output.6.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/tarid/tcf21_rs2329429.pmplot.pdf \
../figures/tarid/tcf21_rs2329429.pvalue_vs_sample_size.pdf 


# make PM plot for TCF21 and rs2327429 (sub sample):
grep "ENSG00000118526.6_6_134209837_T_C_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.6.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/tarid/tcf21_rs2329429.size52.pmplot.pdf none


# make PM plot for TMEM120B and rs148608463 (full sample): 
grep "ENSG00000188735.8_12_121413027_G_A_b37" ../processed_data/160805/metasoft_output/metasoft_output.12.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/tarid/tmem120b_rs148608463.pmplot.pdf \
../figures/tarid/tmem120b_rs148608463.pvalue_vs_sample_size.pdf 


# make PM plot for TMEM120B and rs148608463 (subsample sample):
grep "ENSG00000188735.8_12_121413027_G_A_b37" ../processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.12.mcmc.txt | \
Rscript tarid/pm_plot.R ../figures/tarid/tmem120b_rs148608463.size52.pmplot.pdf none


#### end tarid




#----------------- eQTL and ATACseq ----------# 
# Setup: 
mkdir eqtl_and_atacseq ../processed_data/eqtl_and_atacseq


# Calculate the percentage of eQTLs within ATACseq regions stratified by p-value: 
Rscript eqtl_and_atacseq/overlap_eqtl_and_atacseq.R


# Calculate the percentage of eQTLs within ATACseq regions stratified by p-value and HCASMC-specificity score: 
Rscript eqtl_and_atacseq/overlap_hcasmc_specific_eqtl_and_atacseq.svalue.R # s-value model
Rscript eqtl_and_atacseq/overlap_hcasmc_specific_eqtl_and_atacseq.entropy.R # entropy model


# Link eQTL specificity score: 
ln ../processed_data/hcasmc_specific_eqtl/eqtl_specificity_index/* ../processed_data/eqtl_and_atacseq/


# Calculate eQTL and ATACseq overlap stratified by specificity: 
Rscript eqtl_and_atacseq/overlap_hcasmc_specific_eqtl_and_atacseq.entropy.R


# Gregor: 
Rscript eqtl_and_atacseq/gregor/LD_proxy.R
Rscript eqtl_and_atacseq/gregor/overlap_enrichment.R


#----------------- Brian's genes ---------------# 
# Setup: 
mkdir -p ../figures/brian/locuszoom brian ../processed_data/brian

# make locuszoom for AHR: 
locuszoom --metal ../processed_data/rasqual/output/ENSG00000106546.8_AHR.pval.txt --pvalcol pval --markercol rsid --refsnp rs10265174 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (AHR)" --prefix ../figures/brian/locuszoom/AHR_eQTL
locuszoom --metal ../processed_data/rasqual/output/ENSG00000106546.8_AHR.pval.txt --pvalcol pval --markercol rsid --refsnp rs10265174 --flank 500KB ymax=5 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (AHR)" --prefix ../figures/brian/locuszoom/AHR_eQTL

# MMP1:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1799750 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1799750 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1711437 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs751547 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000196611.4_MMP1.pval.txt --pvalcol pval --markercol rsid --refsnp rs1784405 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000196611.4_MMP1)" --prefix ../figures/brian/locuszoom/ENSG00000196611.4_MMP1_eqtl

# CYP1B1:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000138061.7_CYP1B1.pval.txt --pvalcol pval --markercol rsid --refsnp rs4670837 --flank 1MB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000138061.7_CYP1B1)" --prefix ../figures/brian/locuszoom/ENSG00000138061.7_CYP1B1_eqtl
locuszoom --metal ../processed_data/rasqual/output/ENSG00000138061.7_CYP1B1.pval.txt --pvalcol pval --markercol rsid --refsnp rs162558 --flank 5kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000138061.7_CYP1B1)" --prefix ../figures/brian/locuszoom/ENSG00000138061.7_CYP1B1_eqtl

#### end Brian's genes:


#### metaXcan: 
mkdir metaXcan ../processed_data/metaXcan ../figures/metaXcan

#### end metaXcan


#### PIQ
# TCF21:
cd /srv/persistent/bliu2/tools/piq-single/
Rscript pwmmatch.exact.r common.r \
/srv/persistent/bliu2/shared/jaspar/pfm_all.txt 820 \
$processed_data/footprint/piq/motif.matches/

Rscript bam2rdata.r common.r \
$processed_data/footprint/piq/bam/fbs.RData \
$data/atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam \
$data/atacseq/fbs/2305/out/align/rep2/CA2305-FBS_S2_concat_R1_001.trim.PE2SE.bam \
$data/atacseq/fbs/2305/out/align/rep3/CA2305-FBS_S3_concat_R1_001.trim.PE2SE.bam

Rscript pertf.r common.r \
$processed_data/footprint/piq/motif.matches/ \
$processed_data/footprint/piq/tmp/ \
$processed_data/footprint/piq/fbs_footprint/ \
$processed_data/footprint/piq/bam/fbs.RData 820

# JUN:
Rscript pwmmatch.exact.r common.r \
/srv/persistent/bliu2/shared/jaspar/pfm_all.txt 478 \
$processed_data/footprint/piq/motif.matches/

Rscript pertf.r common.r \
$processed_data/footprint/piq/motif.matches/ \
$processed_data/footprint/piq/tmp/ \
$processed_data/footprint/piq/fbs_footprint/ \
$processed_data/footprint/piq/bam/fbs.RData 478


#### Metabochip eQTL lookup: 
mkdir metabochip ../processed_data/metabochip
parallel "grep '{}' ../processed_data/rasqual/output/*pval.txt | cut -d: -f2 > ../processed_data/metabochip/{}.txt" :::: <(cut -f1 ../data/gwas/metabochip_novel_lead_variants.txt | grep -v markername)
cat ../processed_data/metabochip/*.txt | awk 'BEGIN{OFS="\t"}{if (exp(log(10)*$10)<0.05) print $1,$2,$3,$4,$26,exp(log(10)*$10)}'


#### fine-mapping:
# REST:
locuszoom --metal ../processed_data/rasqual/output/ENSG00000084093.11_REST.pval.txt --pvalcol pval --markercol rsid --refsnp rs17087335 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000084093.11_REST)" --no-date --prefix ../figures/locuszoom/ENSG00000084093.11_REST_eqtl_chr4:57700000-57900000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs17087335 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000084093.11_REST)" --no-date --prefix ../figures/locuszoom/ENSG00000084093.11_REST_gwas_chr4:57700000-57900000

locuszoom --metal ../processed_data/rasqual/output/ENSG00000084093.11_REST.pval.txt --pvalcol pval --markercol rsid --refsnp rs66790703 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_eqtl_chr4:57700000-57900000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs66790703 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_gwas_chr4:57700000-57900000

locuszoom --metal ../processed_data/rasqual/output/ENSG00000084093.11_REST.pval.txt --pvalcol pval --markercol rsid --refsnp rs56155140 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_eqtl_chr4:57700000-57900000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs56155140 --chr chr4 --start 57700000 --end 57900000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000084093.11_REST)" --no-date --plotonly --prefix ../figures/locuszoom/ENSG00000084093.11_REST_gwas_chr4:57700000-57900000

# MAP3K7CL: 
locuszoom --metal ../processed_data/rasqual/output/ENSG00000156265.11_MAP3K7CL.pval.txt --pvalcol pval --markercol rsid --refsnp rs11911017 --chr chr21 --start 30500000 --end 30700000 --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000156265.11_MAP3K7CL)" --no-date --prefix ../figures/locuszoom/ENSG00000156265.11_MAP3K7CL_eqtl_chr21:30500000-30700000
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs11911017 --chr chr21 --start 30500000 --end 30700000 --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000156265.11_MAP3K7CL)" --no-date --prefix ../figures/locuszoom/ENSG00000156265.11_MAP3K7CL_gwas_chr21:30500000-30700000

# COL4A2: 
locuszoom --metal ../processed_data/rasqual/output/ENSG00000134871.13_COL4A2.pval.txt --pvalcol pval --markercol rsid --refsnp rs11838776 --flank 500kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_eqtl
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs11838776 --flank 500kB --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_gwas

locuszoom --metal ../processed_data/rasqual/output/ENSG00000134871.13_COL4A2.pval.txt --pvalcol pval --markercol rsid --refsnp rs11838776 --flank 50kB --source 1000G_March2012 --build hg19 --pop EUR title="eQTL (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_eqtl
locuszoom --metal $data/gwas/cad.add.160614.website.metal.txt --pvalcol p_dgc --markercol markername --refsnp rs11838776 --flank 50kB --source 1000G_March2012 --build hg19 --pop EUR  title="GWAS (ENSG00000134871.13_COL4A2)" --no-date --prefix ../figures/locuszoom/ENSG00000134871.13_COL4A2_gwas


#----------------- MPRA array ----------------# 
# setup: 
mkdir mpra ../processed_data/mpra


# find LD SNPs with rAggr: 
cat ../data/gwas/gwas.cad.mi.recessive.metabochip.txt | sort -n -k2 -n -k3 | uniq | awk '{print "chr"$2":"$3}' > ../processed_data/mpra/rAggr/rAggrInput.txt
# On http://raggr.usc.edu/, select 1000 Genomes Phase 3, East Asian, European, and South Asian, MAF = 0.01, r2 range [0.6,1], Max Distance 500kb, other parameters are left as default. 


# Remove duplicated variants in rAggr output, and select variants with r2 > 0.8 (with the lead variants): 
Rscript mpra/rAggr2uniq_rAggr.R ../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.6_1kgp3.csv ../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.8_1kgp3.uniq.txt 0.8


# select tested genes for each GWAS proxy SNPs: 
mkdir ../processed_data/mpra/naive_overlap
parallel --xapply "grep -P 'chr{1}' ../processed_data/rasqual/output/*pval.txt | cut -d: -f2 > ../processed_data/mpra/naive_overlap/{2}.txt" :::: <(cat ../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.8_1kgp3.uniq.txt | grep -v markername | cut -f2,3) :::: <(cat ../processed_data/mpra/rAggr/maf0.01_dist500kb_rsq0.8_1kgp3.uniq.txt | grep -v markername | cut -f1)


# select eQTLs for all GWAS proxy SNPs: 
mkdir ../processed_data/mpra/naive_overlap_eQTL/
cat ../processed_data/mpra/naive_overlap/*.txt | awk 'BEGIN{OFS="\t";print "fid","rsid","chr","pos","ref","alt","af","pi","n_rsnps","rsq_rsnp","pval","padj","rank"}{if (exp(log(10)*$10)<0.05) print $1,$2,$3,$4,$5,$6,$7,$12,$18,$25,$26,exp(log(10)*$10),$27}' | sort -g -k10 > ../processed_data/mpra/naive_overlap_eQTL/gwas_eQTL_naive_overlap.txt


# QC on colocalized eQTLs: 
Rscript mpra/qc_naive_overlap.R ../processed_data/mpra/naive_overlap_eQTL/gwas_eQTL_naive_overlap.txt ../figures/mpra/naive_overlap/


# Select colocalized genes and compute the number of putative variants for each gene:
mkdir ../processed_data/mpra/eCAVIAR
bash mpra/select_eCAVIAR_genes.sh ../processed_data/eCAVIAR/eCAVIAR_output_rasqual_111616/ ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_genes.txt


# Plot the distribution of putative eQTL/GWAS/intersection/union variants in the 95% confidence set for colocalized/null genes:
mkdir ../figures/mpra/eCAVIAR
bash mpra/qc_eCAVIAR_genes.R ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_genes.txt ../figures/mpra/eCAVIAR/qc_eCAVIAR_genes.pdf


# Select colocalized variants selected by eCAVIAR:
bash mpra/select_eCAVIAR_variants.sh ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_genes.txt ../processed_data/eCAVIAR/eCAVIAR_output_rasqual_111616/ /srv/scratch/bliu2/concat_eCAVIAR_variants/ ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.bed


# Add RASQUAL statistics, conservation, TFBS, open chromatin, CADD score to putative variants selected by eCAVIAR: 
bash mpra/add_features.sh ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.bed ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.allFeatures.bed


# Plot the distribution of annotations (TFBS, CADD score, etc) for variants selected by eCAVIAR:
Rscript mpra/qc_eCAVIAR_variants.R ../processed_data/mpra/eCAVIAR/eCAVIAR_colocalized_variants.allFeatures.bed ../figures/mpra/eCAVIAR/qc_eCAVIAR_variants.pdf


# Subset to columns to provide to Nathan Abell: 
Rscript mpra/format_eCAVIAR_variants.R
Rscript mpra/format_naive_overlap.R


# Choose cell type for MPRA:
Rscript mpra/choose_cell_type.R
Rscript mpra/choose_cell_type.expanded.R
 

#------------------ Feature construction ----------------
# Convert VCF into BEDs:
mkdir ../data/variantBeds feature_construction
for i in {1..22}; do zcat ../data/joint3/phased_imputed/phased_and_imputed.chr$i.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz | grep -v '#' | awk 'BEGIN{FS="\t|;"; OFS="\t"}{print $1,$2-1,$2,$3,$4,$5,$10}' | sed 's/AF=//' | sort -k2,2n > ../data/variantBeds/chr$i.bed; done


# Add CADD score to variant BEDs (output will be gzipped):
mkdir ../data/features 
parallel -j10 bash feature_construction/add.cadd.scores.sh ../data/variantBeds/chr{}.bed ../data/features/chr{}.withCADD.bed ::: {1..22}


# process HCASMC ATACseq file for 2305:
cut -f1,2,3,7 ../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.narrowPeak | sort -k1,1 -k2,2n | gzip > /srv/scratch/bliu2/2305_ppr.IDR0.1.filt.narrowPeak.gz


# Add ATACseq, ENCODE, 1000 Genomes, CADD TF features to variants BEDs: 
parallel -j10 bash feature_construction/add.features.sh ../data/features/chr{}.withCADD ../data/features/chr{}.allFeatures.bed.gz &> ../logs/add.features.log ::: {1..22}


#----------------- Chromatin annotation ------------------
# Setup: 
mkdir -p chromatin ../processed_data/chromatin ../figures/chromatin ../data/chromatin/fastq


# Prepare chromatin peaks files: 
bash chromatin/valk2durga.sh 


# Install ChIPseq pipeline from Kundaje lab and ChromHMM: 
bash chromatin/install.sh 


# Run ChIPseq pipeline on H3K4me1, H3K4me3, H3K27me3, and H3K27ac data: 
cp chromatin/process_histone_mods.sh /users/bliu2/
bash /users/bliu2/process_histone_mods.sh # run on Nandi!


# Install ChromHMM:



#------------------ Shared ----------------
# Setup: 
mkdir -p shared

# File to map each tissue to a color: 
Rscript shared/tissue_color.R   # shared/tissue_color.txt



#------------------ HCASMC-specific open chromatin ---------
# Setup:
mkdir -p ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/ ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc
ln /srv/persistent/bliu2/HCASMC_eQTL/processed_data/mpra/DHS_expanded/* ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc
ln ../data/atacseq/fbs/*/out/peak/macs2/overlap/*.naive_overlap.narrowPeak.gz ../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc


# Preprocess encode data: 
Rscript hcasmc_specific_open_chromatin/process_encode_data.R
Rscript hcasmc_specific_open_chromatin/process_encode_data_without_taking_max.R

# Find shared peaks: 
python hcasmc_specific_open_chromatin/find_shared_peaks.py


# Get signal at shared peaks:
python hcasmc_specific_open_chromatin/get_signal_at_shared_peaks.py
Rscript hcasmc_specific_open_chromatin/merge_signals.R


# Quantile normalization: 
Rscript hcasmc_specific_open_chromatin/quantile_normalization.R


# Define HCASMC-specific peak score:
Rscript hcasmc_specific_open_chromatin/calc_peak_specificity_index.R


# Intersect specificity with sample peaks: 
Rscript hcasmc_specific_open_chromatin/intersect_specificity_and_peaks.R


# Take the minimum specificity index of each peak: 
Rscript hcasmc_specific_open_chromatin/get_min_index.R


# Assign peak specificity score to raw peaks:
Rscript hcasmc_specific_open_chromatin/assign_raw_peak_specificity.R

# Plot the distribution of number of tissues sharing each peak in HCASMC:
Rscript hcasmc_specific_open_chromatin/plot_distribution_of_tissue_sharing.R


# Perform multi-dimensional scaling: 
Rscript hcasmc_specific_open_chromatin/mds_open_chromatin.R 


# Classifying peaks based on TSS, gene body, intergenic:
Rscript hcasmc_specific_open_chromatin/classify_peaks.R


#--------------- Finemapping ------------
# Setup:
mkdir -p finemap/{enloc,finemap} 

# Finemap: 
bash finemap/finemap/howson_preprocess.sh
Rscript finemap/finemap/howson_data.R
bash finemap/finemap/finemap.sh

Rscript finemap/finemap/rasqual_data.R
Rscript finemap/finemap/howson_rasqual_coloc.R

Rscript finemap/finemap/nikpay_rasqual_coloc.R

Rscript finemap/finemap/atacseq_overlap.R
Rscript finemap/finemap/atacseq_overlap.gviz.R

Rscript finemap/finemap/motif_match.R

Rscript finemap/finemap/tarid_tcf21_coloc.R


#--------------- GWAS gene overlap ------------# 
Rscript gwas_gene_overlap/ldscore_regression/tissue_specific_gene.R
Rscript gwas_gene_overlap/ldscore_regression/tissue_specific_snp_annotation.R
bash gwas_gene_overlap/ldscore_regression/gwas_sumstats.sh

bash gwas_gene_overlap/ldscore_regression/ldscore.sh
bash gwas_gene_overlap/ldscore_regression/ldscore_merged.sh

bash gwas_gene_overlap/ldscore_regression/partition_heritability.sh
bash gwas_gene_overlap/ldscore_regression/partition_heritability_merged.sh

