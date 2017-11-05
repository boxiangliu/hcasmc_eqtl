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


#### 160724:
# obj: re-map sQTL using WASP files:
# setup: 
mkdir $scripts/160724 $figures/160724 $processed_data/160724

# run leafcutter: 
cp $scripts/160627/run_bam2junc.sh $scripts/160724/

bash $scripts/160724/run_bam2junc.sh

