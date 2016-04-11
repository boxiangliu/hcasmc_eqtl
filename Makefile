# HRC phase checkup: 
dummy: 
	# Checking your VCF is valid
	bcftools view geno.vcf.gz >/dev/null
	# Making sure your VCF is indexable by bcftools
	bcftools index geno.vcf.gz
	# Keeping only biallelic SNPs and Indels
	bcftools view -m2 -M2 -v snps,indels geno.vcf.gz -Oz -o out.vcf.gz
	mv out.vcf.gz geno.vcf.gz
	# Making sure your VCF contains only genotypes
	bcftools annotate -x FORMAT geno.vcf.gz -Oz -o out.vcf.gz
	mv out.vcf.gz geno.vcf.gz
	# Making sure your VCF does not contain any missing genotypes
	bcftools view -g ^miss geno.vcf.gz -Oz -o out.vcf.gz
	mv out.vcf.gz geno.vcf.gz
	# Subsetting your VCF down to the reference panel site list
	bcftools index geno.vcf.gz
	bcftools index var_list.vcf.gz
	bcftools isec -n=2 -w1 geno.vcf.gz var_list.vcf.gz -Oz -o out.vcf.gz
	mv out.vcf.gz geno.vcf.gz
	# Making sure your VCF contains the chromosomes in the variant list
	bcftools index geno.vcf.gz
	for i in $(seq 1 22); do 
		bcftools view -G -r $i geno.vcf.gz >/dev/null
	done

# imputation: 
input19 = ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.vcf.gzary
bed19 = ../processed_data/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr20.bed
freq19 = ../processed_data/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr20.freq
figure19 = ../figures/plot_percentage_in_1000G.pdf
reference_hap19 = $(SHARED)/haplotype_reference/1000G_phase3/1000GP_Phase3/1000GP_Phase3_chr20.hap.gz 
reference_legend19 = $(SHARED)/haplotype_reference/1000G_phase3/1000GP_Phase3/1000GP_Phase3_chr20.legend.gz
reference_sample19 = $(SHARED)/haplotype_reference/1000G_phase3/1000GP_Phase3/1000GP_Phase3.sample
genetic_map19 = $(SHARED)/haplotype_reference/1000G_phase3/1000GP_Phase3/genetic_map_chr20_combined_b37.txt
log19 = ../processed_data/imputation/alignments.snp.strand
phased19 = ../processed_data/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr20.phased.haps
imputed19 = ../processed_data/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr20.phased.imputed
# create plink bed file: 
$(bed19): $(input19)
	if [ ! -d $(dir $(bed19)) ]; then mkdir $(dir $(bed19)); fi
	plink \
	--vcf $(input19) \
	--keep-allele-order \
	--make-bed \
	--freq \
	--chr 20 \
	--out $(bed19:.bed=)


# check strand alignment:
$(log19): $(bed19)
	shapeit -check --input-bed $(bed19:.bed=) --input-ref $(reference_hap19) $(reference_legend19) $(reference_sample19) --output-log $(log19:.snp.strand=) > log/check_strand_alignment.log
	less ../processed_data/imputation/alignments.snp.strand


# plot percentage in 1000G vs freq in main panel:
$(figure19): $(freq19) $(log19)
	Rscript plot_percentage_in_1000G.R $(freq19) $(log19) $(figure19)



# pre-phasing: 
$(phased19):$(log19)
	shapeit -B $(bed19:.bed=) -M $(genetic_map19) -R $(reference_hap19) $(reference_legend19) $(reference_sample19) -O $(phased19:.haps=) --exclude-snp $(log19).exclude --thread 12 --window 0.5
	# log files: 
	# log/shapeit_09042016_18h26m07s_1525f352-0e83-4e99-b145-3b758be7e262.{log,ind.mm,snp.mm}


# imputation: 
$(imputed19):$(phased19)
	bash impute.sh $(phased19) $(reference_hap19) $(reference_legend19) $(genetic_map19) $(imputed19) 63025520
	# log file: 
	# ../processed_data/imputation/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr20.phased.imputed.*_summary

# imputation pipeline: 
.imputation_pipeline.sh.done:$(input19)
	bash imputation_pipeline.sh $(input19)
# pre-imputation QC: 
# input8 = ../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
# input8 = ../data/joint/recalibrated_variants.GRCh37.pass.id.vcf.gzary
input8 = ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.vcf.gz
output8 = ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.vcf.gzary
exclude = ../processed_data/discordant_loci_ID.txt
geno_miss_thresh = 0.5 
hwe = 1e-12 
maf = 0.01 # 1 allele in 67 individuals
$(output8): $(input8)
	plink \
	--vcf $(input8) \
	--keep-allele-order \
	--make-bed \
	--maf $(maf) \
	--geno $(geno_miss_thresh) \
	--hwe $(hwe) \
	--recode vcf bgz \
	--out $(output8:.vcf.gzary=) > log/pre_imputation_qc.log

	# --exclude $(exclude) \
	tabix -p vcf $(output8)

# RNA-WGS match: 
input17 = ../data/joint/recalibrated_variants.vcf.gz
variants17 = ../processed_data/rna_wgs_match/variants.chr1_chr11_chr22.bed
mpileup_dir17 = ../processed_data/rna_wgs_match/mpileup
sample_list17 = rna_wgs_match.pileup.sample_list.txt
sample_list17_2 = rna_wgs_match.variant_counts.sample_list.txt
figure17 = ../figures/rna_wgs_match.pdf
$(variants17): $(input17)
	# extract variant sites on chr22:
	zcat $(input17) | awk 'BEGIN {OFS = "\t"} {if ( ($$1 == "chr1" || $$1 == "chr11" || $$1 == "chr22") && $$7 == "PASS") print $$1,$$2-1,$$2}' > $(variants17)

.rna_wgs_match.mpileup.sh.done: $(variants17)
	# samtools mpileup at variant sites: 
	bash rna_wgs_match.mpileup.sh ../data/rnaseq/alignments/ $(sample_list17) $(mpileup_dir17) $(variants17) 
	
.rna_wgs_match.variant_counts.sh.done: .rna_wgs_match.mpileup.sh.done $(sample_list17_2)
	# convert mpileup to variant counts: 
	bash rna_wgs_match.variant_counts.sh ../processed_data/rna_wgs_match/mpileup/ $(sample_list17_2) ../processed_data/rna_wgs_match/variant_count/


split:split
	# extract genotype from WGS VCF file: 
	vcftools --gzvcf $(input17) --extract-FORMAT-info GT --out $(input17:.vcf.gz=)

	# check RNA-WGS concordance: 
	Rscript rna_wgs_match.R $(input17:.vcf.gz=.GT.FORMAT.chr22) ../processed_data/rna_wgs_match/variant_count/ rna_wgs_match.R.sample_list.txt

# copy RNAseq alignments from valk: 
input18 = copy_rnaseq_from_valk.sample_list.txt
output18 = .copy_rnaseq_from_valk.done
$(output18):$(input18)
	bash copy_rnaseq_from_valk.sh /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments/ $(input18)


# compare duplicate samples: 
input15 = ../processed_data/recalibrated_variants.GRCh37.tsv
output15_1 = ../processed_data/discordant_loci_distribution.tsv 
output15_2 = ../processed_data/discordant_loci_ID.txt
$(output15_2): $(input15)
	Rscript compare_duplicate_samples.R $(input15) $(output15_1) $(output15_2)




# missingness rate using plink: 
input14 = ../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
output14 = ../data/joint/recalibrated_variants.GRCh37.pass
figure14 = ../figures/missingness.pdf
$(figure14): $(input14)
	plink --vcf $(input14) --keep-allele-order --missing --out ../data/joint/recalibrated_variants.GRCh37.pass
	Rscript plot_missingness_rate.R ../data/joint/recalibrated_variants.GRCh37.pass ../figures/missingness.pdf



# convert VCF file to table: 
input13 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output13 =  ../processed_data/recalibrated_variants.GRCh37.tsv
$(output13): $(input13)
	java -jar $(GATK) -R $(GRCh37) -T VariantsToTable -V $(input13) --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -GF GT -o ../processed_data/recalibrated_variants.GRCh37.tsv 


# set missing variant ids: 
# TODO: normalized the vcf file first
# input16 = ../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
# output16 = ../data/joint/recalibrated_variants.GRCh37.pass.id.vcf.gzary
# $(output16): $(input16)
# 	plink \
# 	--vcf $(input16) \
# 	--keep-allele-order \
# 	--set-missing-var-ids @_#_\$$1_\$$2 \
# 	--recode vcf bgz \
# 	--out $(output16:vcf.gzary=)
# 	tabix -p vcf $(output16)

# normalize and left-align indels, check if REF allele matches the reference genome: 
input21 = ../data/joint/recalibrated_variants.GRCh37.biallelic.vcf.gz
output21 = ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.vcf.gz
$(output21): $(input21)
	bcftools view -Ou -f .,PASS $(input21) | bcftools norm -Ou -f $(GRCh37) | bcftools annotate -Oz -I +'%CHROM\_%POS\_%REF\_%ALT' -o $(output21)
	# Reference allele mismatch at Y:2649246 .. 'N' vs 'C'
	tabix -p vcf $(output21)


# keep only PASS'ed variants: 
input3 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output3 = ../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
$(output3): $(input3)
	java -jar /usr/bin/GenomeAnalysisTK.jar -R $(GRCh37) -T SelectVariants -V $(input3) -o $(output3:vcf.gz=vcf) --excludeFiltered 2> log/SelectVariants.log 
	bgzip $(output3:vcf.gz=vcf)
	tabix -p vcf $(output3)


# plot ancestry proportions and PCs:  
aims = ../processed_data/top_1000_aims.bed
input7 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.vcf.gz
1000G_phase3v5a = /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
merged = ../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.vcf.gz
1000G_panel = /srv/persistent/bliu2/shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel
figure7_1 = ../figures/ancestry_proportions.pdf
figure7_2 = ../figures/ancestry_PCs.pdf
$(figure7_1): $(input7) 
	# create $(aims)
	tail -n +2 /srv/persistent/bliu2/ancestry/AIMS_selection/AIMs/five_superpopulations/with_heterogeneity_filter/filter_SAS/propmultiIn_0.0/AFR_AMR_EAS_EUR_SAS_with_heterogeneity_filter500k_1000.aims | awk 'BEGIN {OFS = "\t"}{print $$2,$$3-1,$$3}' > $(aims) 
	
	# subset at AIM loci: 
	tabix -h $(input7) -B $(aims) > $(input7:vcf.gz=1000aims.vcf); bgzip $(input7:vcf.gz=1000aims.vcf); tabix -p vcf $(input7:vcf.gz=1000aims.vcf).gz
	tabix -h $(1000G_phase3v5a) -B $(aims) > $(1000G_phase3v5a:vcf.gz=1000aims.vcf); bgzip $(1000G_phase3v5a:vcf.gz=1000aims.vcf); tabix -p vcf $(1000G_phase3v5a:vcf.gz=1000aims.vcf).gz 
	bcftools merge $(input7:vcf.gz=1000aims.vcf.gz) $(1000G_phase3v5a:vcf.gz=1000aims.vcf.gz) -Oz > $(merged)
	vcftools --gzvcf $(merged) --extract-FORMAT-info GT --out $(merged:.gz=)

	# remove multi-allelic and missing loci:
	bcftools view -m2 -M2 -v snps,indels -g ^miss $(merged) -Oz -o $(merged:vcf.gz=biallelic.nomissing.vcf.gz)

	# remove SAS population: 
	awk '{if ($$3 == "SAS") print $$1}' $(1000G_panel) > ../processed_data/SAS_individuals.txt
	bcftools view -S ^../processed_data/SAS_individuals.txt $(merged:vcf.gz=biallelic.nomissing.vcf.gz) -Oz -o $(merged:vcf.gz=biallelic.nomissing.noSAS.vcf.gz)
	
	# calculate using ADMIXTURE: 
	plink --vcf $(merged:vcf.gz=biallelic.nomissing.noSAS.vcf.gz) --keep-allele-order --make-bed --out $(merged:vcf.gz=biallelic.nomissing.noSAS)
	admixture $(merged:vcf.gz=biallelic.nomissing.noSAS.bed) 4
	mv $(notdir $(merged:vcf.gz=biallelic.nomissing.noSAS.4.P)) ../processed_data/
	mv $(notdir $(merged:vcf.gz=biallelic.nomissing.noSAS.4.Q)) ../processed_data/

	# plot ancestry proportion with R: 
	Rscript plot_ancestry_proportions.R ../processed_data/$(notdir $(merged:vcf.gz=biallelic.nomissing.noSAS.4.Q)) ../processed_data/$(notdir $(merged:vcf.gz=biallelic.nomissing.noSAS.fam)) $(figure7_1)

	# plot ancestry PCs with R:
	Rscript plot_ancestry_PCs.R ../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.vcf.GT.FORMAT $(figure7_2)


# visualize WGS coverage: 
input12 = /mnt/data/WGS_HCASMC/sample_list.txt
output12 = .visualize_wgs_coverage.done
$(output12): $(input12)
	bash visualize_wgs_coverage.sh /mnt/data/WGS_HCASMC/ $(input12) $(output12)

# plot genotype correlation: 
input9 = ../data/joint/recalibrated_variants.GRCh37.biallelic.nomissing.vcf.gz
figure9 = ../figures/genotype_correlation_heatmap_full_genotype.pdf 
$(figure9):$(input9)
	# extract genotypes: 
	vcftools --gzvcf $(input9) --extract-FORMAT-info GT --out $(input9:.vcf.gz=)
	# plot genotype correlation: 
	Rscript plot_genotype_correlation.R ./ $(input:.gz=.GT.FORMAT) $(figure9)
	# variant different between 150328 and 59386145:
	tabix -h ../data/joint/recalibrated_variants.GRCh37.vcf.gz -B ../processed_data/variants_different_between_150328_and_59386145.bed > ../processed_data/recalibrated_variants.GRCh37.vcf.variants_different_between_150328_and_59386145

# remove multiallelic loci and missing genotypes: 
input10 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output10 = $(input10:vcf.gz=biallelic.nomissing.vcf.gz)
$(output10): $(input10)
	# Keeping only biallelic SNPs and Indels; removing missing genotypes. 
	bcftools view -m2 -M2 -v snps,indels -g ^miss $(input10) -Oz -o $(output10)
	tabix -p vcf $(output10)

# remove multiallelic loci: 
input20 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output20 = ../data/joint/recalibrated_variants.GRCh37.biallelic.vcf.gz
$(output20): $(input20)
	# Keeping only biallelic SNPs and Indels 
	bcftools view -m2 -M2 -v snps,indels $(input20) -Oz -o $(output20)
	tabix -p vcf $(output20)

	# what percentage are multi-allelic sites: 
	bcftools view -H $(input20) | wc -l  # 26604334
	bcftools view -H $(output20) | wc -l # 25275762
	# a=26604334
	# b=25275762
	# (a-b)/a = 0.04993818

# plot plink/seq istats files:
input6 = ../processed_data/recalibrated_variants.vcf.pseq.istats
figure6 = ../figures/pseq_istats_recalibrated_variants.pdf
$(figure6): $(input6)
	Rscript plot_pseq_stats.R ../processed_data

# vcf quality control using plink/seq: 
../processed_data/recalibrated_variants.vcf.pseq.vstats: /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf
	pseq /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf v-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.vcf.pseq.vstats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf i-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.vcf.pseq.istats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf v-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.vstats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf i-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.istats


# vcf quality control using QC3: 
input4 = /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf
output4 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.QC3
$(output4): $(input4)
	perl /srv/persistent/bliu2/tools/QC3/qc3.pl -m v -i $(input4) -a /srv/persistent/bliu2/tools/annovar/humandb -o $(output4)

input5 = /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf
output5 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.QC3
$(output5): $(input5)
	perl /srv/persistent/bliu2/tools/QC3/qc3.pl -m v -i $(input5) -a /srv/persistent/bliu2/tools/annovar/humandb -o $(output5)

# genotype refinement: 
/mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.gz: /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.gz
	python genotype_refinement.py /mnt/data/WGS_HCASMC/joint/ recalibrated_variants.GRCh37.vcf.gz 2>> log/genotype_refinement.log
	cd /mnt/data/WGS_HCASMC/joint/; mkdir genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.idx genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.vcf.idx genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.vcf genotype_refinement_temp

# create dict for GRCh37 fasta: 
input2 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa 
output2 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.dict
$(output2): $(input2) 
	java -jar /software/picard-tools/1.92/CreateSequenceDictionary.jar R=$(input2) O=$(output2)

# convert vcf from hg19 to GRCh37: 
input1 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.gz
output1 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.vcf
$(output1): $(input1)
	zcat $(input1) | sed -e "s/^chr//" -e "s/contig=<ID=chr/contig=<ID=/" > $(output1)
	bgzip $(output1)
	tabix -p vcf $(output1:vcf=vcf.gz)

# concatenate 1000G vcfs: 
/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz: 
	python concat_1000G_vcfs.py /srv/persistent/bliu2/shared/1000genomes/phase3v5a sample_list.txt ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf

# copy BAM and VCF to valk: 
transfer_to_valk.py.done: 
	python transfer_to_valk.py /mnt/data/WGS_HCASMC sample_list.txt

.copy_bam_and_vcf_to_valk.sh.done: 
	bash copy_bam_and_vcf_to_valk.sh copy_bam_and_vcf_to_valk.sample_list.txt

# joint genotype calling: 
/mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf:
	python genotype_gvcfs.py /mnt/data/WGS_HCASMC/joint/ ../sample_list.txt 2> log/genotype_gvcfs.20160319.log 
	# relocate temporary files:
	cd /mnt/data/WGS_HCASMC/joint
	mkdir genotype_gvcfs_temp
	mv recalibrated_snps_raw_indels.vcf* genotype_gvcfs_temp
	mv recalibrated_variants.vcf* genotype_gvcfs_temp
	mv recalibrate_INDEL* genotype_gvcfs_temp
	mv recalibrate_SNP* genotype_gvcfs_temp

# generate md5sum: 
output11 = /mnt/data/WGS_HCASMC/md5sum
$(output11):
	cd /mnt/data/WGS_HCASMC/
	md5sum */*{fastq.gz,vcf.gz,vcf.gz.tbi,bai,bam} > md5sum

# run 1369: 
# /mnt/data/WGS_HCASMC/1369/raw_variants.g.vcf.gz: 
# 	bash wgs_pipeline/wgs_pipeline.1369.160312.sh


# check finished samples: 
/mnt/data/WGS_HCASMC/check_finished_sample.py.done:
	# python check_finished_sample.py /mnt/data/WGS_HCASMC sample_list_35samples.txt
	python check_finished_sample.py /mnt/data/WGS_HCASMC sample_list.txt

# compress g.vcf files: 
/mnt/data/WGS_HCASMC/compress_gVCF.py.done: 
	python compress_gVCF.py /mnt/data/WGS_HCASMC sample_list_3samples.txt

HG19 = /srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
GATK =  /usr/bin/GenomeAnalysisTK.jar
GRCh37 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa
SHARED = /srv/persistent/bliu2/shared/