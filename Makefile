HG19 = /srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
GATK =  /usr/bin/GenomeAnalysisTK.jar
GRCh37 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa
SHARED = /srv/persistent/bliu2/shared/


# compress g.vcf files: 
/mnt/data/WGS_HCASMC/compress_gVCF.py.done: 
	python compress_gVCF.py /mnt/data/WGS_HCASMC sample_list_3samples.txt


# check finished samples: 
/mnt/data/WGS_HCASMC/check_finished_sample.py.done:
	# python check_finished_sample.py /mnt/data/WGS_HCASMC sample_list_35samples.txt
	python check_finished_sample.py /mnt/data/WGS_HCASMC sample_list.txt


# run 1369: 
# /mnt/data/WGS_HCASMC/1369/raw_variants.g.vcf.gz: 
# 	bash wgs_pipeline/wgs_pipeline.1369.160312.sh


# generate md5sum: 
output11 = /mnt/data/WGS_HCASMC/md5sum
$(output11):
	cd /mnt/data/WGS_HCASMC/
	md5sum */*{fastq.gz,vcf.gz,vcf.gz.tbi,bai,bam} > md5sum


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


#----- copy BAM and VCF to valk -----# 
# transfer_to_valk.py.done: 
# 	python transfer_to_valk.py /mnt/data/WGS_HCASMC sample_list.txt


.copy_bam_and_vcf_to_valk.sh.done: 
	bash copy_bam_and_vcf_to_valk.sh copy_bam_and_vcf_to_valk.sample_list.txt


# concatenate 1000G vcfs: 
/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz: 
	python concat_1000G_vcfs.py /srv/persistent/bliu2/shared/1000genomes/phase3v5a sample_list.txt ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf


# convert vcf from hg19 to GRCh37: 
input1 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.gz
output1 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.vcf
$(output1): $(input1)
	zcat $(input1) | sed -e "s/^chr//" -e "s/contig=<ID=chr/contig=<ID=/" > $(output1)
	bgzip $(output1)
	tabix -p vcf $(output1:vcf=vcf.gz)


# create dict for GRCh37 fasta: 
input2 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa 
output2 = /srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.dict
$(output2): $(input2) 
	java -jar /software/picard-tools/1.92/CreateSequenceDictionary.jar R=$(input2) O=$(output2)


# genotype refinement: 
/mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.gz: /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.gz
	python genotype_refinement.py /mnt/data/WGS_HCASMC/joint/ recalibrated_variants.GRCh37.vcf.gz 2>> log/genotype_refinement.log
	cd /mnt/data/WGS_HCASMC/joint/; mkdir genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.idx genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.vcf.idx genotype_refinement_temp
	mv recalibrated_variants.GRCh37.postCGP.vcf genotype_refinement_temp


# vcf quality control using QC3: 
input4 = /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf
output4 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.QC3
$(output4): $(input4)
	perl /srv/persistent/bliu2/tools/QC3/qc3.pl -m v -i $(input4) -a /srv/persistent/bliu2/tools/annovar/humandb -o $(output4)

input5 = /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf
output5 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.vcf.QC3
$(output5): $(input5)
	perl /srv/persistent/bliu2/tools/QC3/qc3.pl -m v -i $(input5) -a /srv/persistent/bliu2/tools/annovar/humandb -o $(output5)


# vcf quality control using plink/seq: 
../processed_data/recalibrated_variants.vcf.pseq.vstats: /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf
	pseq /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf v-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.vcf.pseq.vstats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_gvcfs_temp/recalibrated_variants.vcf i-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.vcf.pseq.istats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf v-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.vstats
	pseq /mnt/data/WGS_HCASMC/joint/genotype_refinement_temp/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf i-stats --mask reg=chr22 any.filter.ex > ../processed_data/recalibrated_variants.GRCh37.postCGP.GQfiltered.vcf.pseq.istats



# plot plink/seq istats files:
input6 = ../processed_data/recalibrated_variants.vcf.pseq.istats
figure6 = ../figures/pseq_istats_recalibrated_variants.pdf
$(figure6): $(input6)
	Rscript plot_pseq_stats.R ../processed_data


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


# remove multiallelic loci and missing genotypes: 
input10 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output10 = $(input10:vcf.gz=biallelic.nomissing.vcf.gz)
$(output10): $(input10)
	# Keeping only biallelic SNPs and Indels; removing missing genotypes. 
	bcftools view -m2 -M2 -v snps,indels -g ^miss $(input10) -Oz -o $(output10)
	tabix -p vcf $(output10)


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


# visualize WGS coverage: 
input12 = /mnt/data/WGS_HCASMC/sample_list.txt
output12 = .visualize_wgs_coverage.done
$(output12): $(input12)
	bash visualize_wgs_coverage.sh /mnt/data/WGS_HCASMC/ $(input12) $(output12)


# plot ancestry proportions and PCs:  
aims = ../processed_data/top_1000_aims.bed
input7 = /mnt/data/WGS_HCASMC/joint/recalibrated_variants.GRCh37.vcf.gz
1000G_phase3v5a = /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
merged = ../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.vcf.gz
1000G_panel = /srv/persistent/bliu2/shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel
figure7_1 = ../figures/ancestry_proportions.pdf
figure7_2 = ../figures/ancestry_PCs.pdf
table7 = ../processed_data/ancestry_proportions.tsv
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
	Rscript plot_ancestry_proportions.R ../processed_data/$(notdir $(merged:vcf.gz=biallelic.nomissing.noSAS.4.Q)) ../processed_data/$(notdir $(merged:vcf.gz=biallelic.nomissing.noSAS.fam)) $(figure7_1) $(table7)

	# plot ancestry PCs with R:
	Rscript plot_ancestry_PCs.R ../processed_data/recalibrated_variants.1000G_phase3v5a.merged.1000aims.vcf.GT.FORMAT $(figure7_2)


# keep only PASS'ed variants: 
input3 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output3 = ../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
$(output3): $(input3)
	java -jar /usr/bin/GenomeAnalysisTK.jar -R $(GRCh37) -T SelectVariants -V $(input3) -o $(output3:vcf.gz=vcf) --excludeFiltered 2> log/SelectVariants.log 
	bgzip $(output3:vcf.gz=vcf)
	tabix -p vcf $(output3)


# normalize and left-align indels, check if REF allele matches the reference genome: 
input21 = ../data/joint/recalibrated_variants.GRCh37.biallelic.vcf.gz
output21 = ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.vcf.gz
$(output21): $(input21)
	bcftools view -Ou -f .,PASS $(input21) | bcftools norm -Ou -f $(GRCh37) | bcftools annotate -Oz -I +'%CHROM\_%POS\_%REF\_%ALT' -o $(output21)
	# Reference allele mismatch at Y:2649246 .. 'N' vs 'C'
	tabix -p vcf $(output21)


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


# convert VCF file to table: 
input13 = ../data/joint/recalibrated_variants.GRCh37.vcf.gz
output13 =  ../processed_data/recalibrated_variants.GRCh37.tsv
$(output13): $(input13)
	java -jar $(GATK) -R $(GRCh37) -T VariantsToTable -V $(input13) --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -GF GT -o ../processed_data/recalibrated_variants.GRCh37.tsv 


# missingness rate using plink: 
input14 = ../data/joint/recalibrated_variants.GRCh37.pass.vcf.gz
output14 = ../data/joint/recalibrated_variants.GRCh37.pass
figure14 = ../figures/missingness.pdf
$(figure14): $(input14)
	plink --vcf $(input14) --keep-allele-order --missing --out ../data/joint/recalibrated_variants.GRCh37.pass
	Rscript plot_missingness_rate.R ../data/joint/recalibrated_variants.GRCh37.pass ../figures/missingness.pdf


# compare duplicate samples: 
input15 = ../processed_data/recalibrated_variants.GRCh37.tsv
output15_1 = ../processed_data/discordant_loci_distribution.tsv 
output15_2 = ../processed_data/discordant_loci_ID.txt
$(output15_2): $(input15)
	Rscript compare_duplicate_samples.R $(input15) $(output15_1) $(output15_2)



#------- copy RNAseq alignments from valk ------# 
input18 = copy_rnaseq_from_valk.sample_list.txt
output18 = .copy_rnaseq_from_valk.done
input18_2 = copy_rnaseq_from_valk_2.sample_list.txt
output18_2 = .copy_rnaseq_from_valk_2.done
$(output18):$(input18)
	bash copy_rnaseq_from_valk.sh /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments/ $(input18)


$(output18_2):$(input18_2)
	bash copy_rnaseq_from_valk_2.sh /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments/ $(input18_2)
	
	cd ../data/rnaseq/alignments
	mkdir 59885590_26425_CGATGT_2
	mv 59885590_26425_CGATGT/{Aligned.out.bam.sort.bam,Aligned.out.bam.sort.bam.bai,Log.final.out,Log.progress.out,Log.out} 59885590_26425_CGATGT_2
	
	mv 59885590_26425_CGATGT_2/Aligned.out.bam.sort.bam 59885590_26425_CGATGT_2/Aligned.out.bam
	mv 59885590_26425_CGATGT_2/Aligned.out.bam.sort.bam.bai 59885590_26425_CGATGT_2/Aligned.out.bam.bai
	
	mv 9071501.8_26436_ATCACG/Aligned.out.bam.sort.bam 9071501.8_26436_ATCACG/Aligned.out.bam
	mv 9071501.8_26436_ATCACG/Aligned.out.bam.sort.bam.bai 9071501.8_26436_ATCACG/Aligned.out.bam.bai
	
	mv 8072501_26426_TGACCA/Aligned.out.bam.sort.bam 8072501_26426_TGACCA/Aligned.out.bam
	mv 8072501_26426_TGACCA/Aligned.out.bam.sort.bam.bai 8072501_26426_TGACCA/Aligned.out.bam.bai

	rsync -avzh --progress bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/HCASMC/Raw_Data/Penn_2nd_round_RNAseq_GOOD_cvrg_8-11-15/2913_26442_GCCAAT/Pass2/{Log.out,Log.final.out,Log.progress.out,Aligned.out.bam.sort.bam,Aligned.out.bam.sort.bam.bai} /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq/alignments/2913_26442_GCCAAT
	mv 2913_26442_GCCAAT/Aligned.out.bam.sort.bam 2913_26442_GCCAAT/Aligned.out.bam
	mv 2913_26442_GCCAAT/Aligned.out.bam.sort.bam.bai 2913_26442_GCCAAT/Aligned.out.bam.bai


#------- link Olga's samples: -------# 
output23 = .link_olgas_samples.sh.done
$(output23):
	bash link_olgas_samples.sh 


#---- check all RNAseq has been mapped ------# 
$(output22):$(input22)
	# on valk: 
	mkdir ../processed_data/check_all_rnaseq_has_mapped/
	scp bosh@valkyr.stanford.edu:/home/mpjanic/Stanford_3rd_round_seqs_GOOD_pools1and2only/AAA-StudyInfo.txt ../processed_data/check_all_rnaseq_has_mapped/Stanford_3rd_round_AAA-StudyInfo.txt
	Rscript check_all_rnaseq_has_mapped.R 
	# sample 2913 is in /home/diskstation/RNAseq/HCASMC/Raw_Data/Penn_2nd_round_RNAseq_GOOD_cvrg_8-11-15/


#------- RNA-WGS match -------#
input17 = ../data/joint/recalibrated_variants.vcf.gz
genotype17 = ../data/joint/recalibrated_variants.GT.FORMAT
variants17 = ../processed_data/rna_wgs_match/variants.chr1_chr11_chr22.bed
mpileup_dir17 = ../processed_data/rna_wgs_match/mpileup
count_dir17= ../processed_data/rna_wgs_match/variant_count/
# sample_list17 = rna_wgs_match.pileup.sample_list.txt
# sample_list17 = rna_wgs_match.pileup.sample_list.2.txt
# sample_list17 = rna_wgs_match.pileup.sample_list.3.txt
sample_list17 = rna_wgs_match.pileup.sample_list.4.txt
# sample_list17_2 = rna_wgs_match.variant_counts.sample_list.txt
# sample_list17_2 = rna_wgs_match.variant_counts.sample_list.2.txt
# sample_list17_2 = rna_wgs_match.variant_counts.sample_list.3.txt
sample_list17_2 = rna_wgs_match.variant_counts.sample_list.4.txt
sample_list17_3 = rna_wgs_match.R.sample_list.txt
figure17 = ../figures/rna_wgs_match.pdf
table17 = ../processed_data/rna_wgs_match.tsv
$(variants17): $(input17)
	# extract variant sites on chr22:
	zcat $(input17) | awk 'BEGIN {OFS = "\t"} {if ( ($$1 == "chr1" || $$1 == "chr11" || $$1 == "chr22") && $$7 == "PASS") print $$1,$$2-1,$$2}' > $(variants17)

.rna_wgs_match.mpileup.sh.done: $(variants17) $(sample_list17)
	# samtools mpileup at variant sites: 
	bash rna_wgs_match.mpileup.sh ../data/rnaseq/alignments/ $(sample_list17) $(mpileup_dir17) $(variants17) 

.rna_wgs_match.variant_counts.sh.done: .rna_wgs_match.mpileup.sh.done $(sample_list17_2)
	# convert mpileup to variant counts: 
	bash rna_wgs_match.variant_counts.sh ../processed_data/rna_wgs_match/mpileup/ $(sample_list17_2) $(count_dir17)


$(figure17): .rna_wgs_match.variant_counts.sh.done $(sample_list17_3)
	# extract genotype from WGS VCF file: 
	# vcftools --gzvcf $(input17) --extract-FORMAT-info GT --out $(genotype17:.GT.FORMAT=)

	# check RNA-WGS concordance: 
	Rscript rna_wgs_match.R $(genotype17) $(count_dir17) $(sample_list17_3) $(figure17) $(table17)



#----- pre-imputation QC ------# 
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

	tabix -p vcf $(output8)


#---- imputation ------# 
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


# re-run failed jobs 
.rerun_failed_imputation_jobs.py.done: 
	python rerun_failed_imputation_jobs.py ../data/joint/
	touch .rerun_failed_imputation_jobs.py.done

	# move files into a directory: 
	mkdir ../data/joint/imputation
	mv ../data/joint/*.imputed.* ../data/joint/imputation
	mv ../data/joint/*chr*.{bed,bim,fam,frq,log,nosex,phased,snp}* ../data/joint/imputation

# run on scg3: 
wd=/srv/persistent/bliu2/tools/
copy_imputation_pipeline_to_scg:
	# scp imputation_pipeline.sh subset_to_chromosome.sh check_strand_alignment.sh plot_percentage_in_1000G.R prephasing.sh impute.sh utils.R bliu2@scg3.stanford.edu:/srv/gsfs0/projects/montgomery/bliu2/HCASMC_eQTL/scripts/
	# scp ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.vcf.gzary bliu2@scg3.stanford.edu:/srv/gsfs0/projects/montgomery/bliu2/HCASMC_eQTL/data/joint/
	# scp -r $(wd)/impute_v2.3.2_x86_64_static $(wd)/plink_1.90_beta3_linux_x86_64 $(wd)/shapeit.v2.r790.RHELS_5.4.dynamic bliu2@scg3.stanford.edu:/srv/gsfs0/projects/montgomery/bliu2/tools/
	# scp -r /srv/persistent/bliu2/shared/haplotype_reference/ bliu2@scg3.stanford.edu:/srv/gsfs0/projects/montgomery/bliu2/shared/haplotype_reference
	# scp utils.R bliu2@scg3.stanford.edu:/srv/gsfs0/projects/montgomery/bliu2/HCASMC_eQTL/scripts/




# concatenate chunks in each chromosome 
# I would like to concatenate the info files at some point too. 
# I wrote code in concat_impute2_output.R (commented out) but it does not work
# as intended because the info files have headers and cannot be simply 
# concatenated.
.concat_impute2_output.py.done:
	Rscript concat_impute2_output.R ../data/joint/imputation/


#----- sex determination -----#
input23= ../data/joint/recalibrated_variants.GT.FORMAT
table23= ../processed_data/sex.tsv
figure23 = ../figures/sex.pdf
$(figure23):$(outpu23)
	Rscript sex_determination.R $(input23) $(figure23) $(table23)


#----- count expression ------#
input_dir24=../data/rnaseq/alignments/
input24 =htseq_count.sample_list.txt
output_dir_prefix24=../data/rnaseq/expression/
.htseq_count.sh.done: $(input24) $(input_dir24)
	bash htseq_count.sh $(input_dir24) $(input24) $(output_dir_prefix24)


#----- detect DNA sample contamination ------#
input25 = ../data/joint/recalibrated_variants.GT.FORMAT
figure25 = ../figures/dna_contamination.pdf
1000G_phase3v5a_chr20 = /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
out_dir25 = ../processed_data/dna_contamination/
eur25 = ../processed_data/dna_contamination/chr20.AF.EUR.2.vcf.gz
output25 = ../processed_data/dna_contamination/1020301
# plot het/hom ratio: 
$(figure25):$(input25) 
	Rscript detect_sample_contamination.R $(input25) $(figure25)

# model based method: 
$(eur25):$(input25)
	# if [ ! -d $(out_dir25) ]; then mkdir -p $(out_dir25); fi
	# bcftools view -Ou -s HG00096 --force-samples --types snps $(1000G_phase3v5a_chr20) | bcftools annotate -Ov -x ^INFO/EUR_AF | sed -e "s/EUR_AF/AF/" -e "s/^20/chr20/" | bcftools view -Oz -o $(eur25)
 	bcftools view -Ov ../processed_data/dna_contamination/chr20.AF.EUR.vcf.gz | sed "s/contig=<ID=/contig=<ID=chr/" | bcftools view -Oz -o ../processed_data/dna_contamination/chr20.AF.EUR.2.vcf.gz
$(output25): $(eur25)
	verifyBamID --vcf $(eur25) --bam /mnt/data/WGS_HCASMC/1020301/recal_reads.bam --chip-none --precise --verbose --minAF 0 --minCallRate 0 --out $(out_dir25)/1020301 


#----- IBD -----# 
input26 = ../data/joint/recalibrated_variants.GT.FORMAT
subset26 = ../data/joint/recalibrated_variants.chr22.GT.FORMAT
figure26 = ../figures/ibd.pdf
$(figure26):$(input26)
	awk '{ if ($$1 == "chr22" || $$1 == "CHROM") print $$0}' $(input26) > $(subset26)


#----- variance stabilization --------
input27 = read_htseq_count.sample_list.txt
table27 = ../processed_data/count_matrix.tsv
$(table27): $(input27)
	Rscript 027_read_htseq_count.R $(input27) $(table27)
	Rscript 027_variance_stabilize.R


#----- 28 prepare genotype ----
# 28.0 post imputation QC: 

# 28.1 convert impute2 output to genotype 
mkdir ../processed_data/028_imputed_genotype
Rscript 028_convert_impute2_output_to_genotype.R


# 28.2 compare imputed genotype with unimputed genotypes: 
zcat ../data/joint/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.vcf.gzary | grep -e ^#CHROM -e ^22 >  ../processed_data/028_imputed_genotype/recalibrated_variants.GRCh37.biallelic.pass.norm.id.hwe.missing.maf.chr22.vcf
Rscript 028_compare_imputed_and_original_genotypes.R


#------ 29 prepare gene expression -----
# 29 set up gene expression for eQTL analysis: 
# during meeting with Trieu and Milos
# we renamed RNAseq samples and associated WGS samples
# create symbolic link to indicate the change in RNAseq sample name:
Rscript 029_link_RNAseq_samples.R


#---- 30 PCA analysis -----
# 30.1 save htseq count files for RNAseq sample in the working set into a tsv file 
cp 027_read_htseq_count.R 030_read_htseq_count.R
# modified 030_read_htseq_count.R
Rscript 030_read_htseq_count.R

# 30.2 compare 9052004_dase and 9052004
Rscript 030_compare_9052004_samples.R 


# 30.3 perform DESeq variance stabilization
cp 027_variance_stabilize.R 030_variance_stabilize.R
# modified 030_variance_stabilize.R
subl 030_variance_stabilize.R
Rscript 030_variance_stabilize.R


# 30.4 make covariates table: 
# don't run:
subl 030_covariates_table.R
Rscript 030_covariates_table.R


# 30.4  perform PCA on RNAseq sample in the working set:
subl 030_expression_PCA.R
Rscript 030_expression_PCA.R


#--- 31 Matrix eQTL ----

# 31.1 getting started with some sample code:
# don't run 
# download example data
mkdir ../processed_data/031_matrix_eqtl_example
wget -P ../processed_data/031_matrix_eqtl_example http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt
wget -P ../processed_data/031_matrix_eqtl_example http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/GE.txt
wget -P ../processed_data/031_matrix_eqtl_example http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt
wget -P ../processed_data/031_matrix_eqtl_example http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/geneloc.txt
wget -P ../processed_data/031_matrix_eqtl_example http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/snpsloc.txt
subl 031_matrix_eqtl_example.R



# 31.2 prepare matrix eQTL genotypes:
# SNP file is already in matrix eQTL format, but 
# their columns were not the same as the working set. So 
# their columns were subsetted to the working set and sorted
# alphanumerically. 
subl 031_prepare_matrix_eQTL_genotype.R
for chr in $(seq 1 22); do
echo $chr
Rscript 031_prepare_matrix_eQTL_genotype.R $chr &
done


# 31.3 prepare matrix eQTL expression: 
# create output directory: 
mkdir ../processed_data/031_prepare_matrix_eQTL_expression
# export DESeq vsd to matrix eQTL expression:
subl 031_prepare_matrix_eQTL_expression.R
Rscript 031_prepare_matrix_eQTL_expression.R 



# 31.4 prepare matrix eQTL covariates 
mkdir ../processed_data/031_prepare_matrix_eQTL_covariate
subl 031_prepare_matrix_eQTL_covariate.R
Rscript 031_prepare_matrix_eQTL_covariate.R

# 31.5 keep snp with MAF > 5%
mkdir ../processed_data/031_subset_genotype_by_maf
subl 031_subset_genotype_by_maf.R
for chr in $(seq 1 22); do
echo $chr
Rscript 031_subset_genotype_by_maf.R $chr &
done


# 31.6 prepare matrix eQTL snp location
mkdir ../processed_data/031_gen_snps_loc
subl 031_gen_snps_loc.R


for chr in $(seq 1 22); do
echo $chr
Rscript 031_gen_snps_loc.R $chr &
done


# 31.7 prepare matrix eQTL gene location 
mkdir ../processed_data/031_gen_gene_loc
subl 031_gen_gene_loc.R
Rscript 031_gen_gene_loc.R
# less ../processed_data/031_gen_gene_loc/gene_loc.txt


# 31.8 find optimal number of PCs:
mkdir ../processed_data/031_find_optimal_PCs
# gather matrix eQTL into one place: 
ln ../processed_data/031_prepare_matrix_eQTL_expression/expression.txt ../processed_data/031_find_optimal_PCs/expression.txt
ln ../processed_data/031_subset_genotype_by_maf/chr20.genotype.txt ../processed_data/031_find_optimal_PCs/chr20.genotype.txt
ln ../processed_data/031_prepare_matrix_eQTL_covariate/covariates.txt ../processed_data/031_find_optimal_PCs/covariates.txt
ln ../processed_data/031_gen_snps_loc/chr20.genotype_loc.txt ../processed_data/031_find_optimal_PCs/chr20.genotype_loc.txt
ln ../processed_data/031_gen_gene_loc/gene_loc.txt ../processed_data/031_find_optimal_PCs/gene_loc.txt
awk '{ if ($2=="20" || $2=="chr") print $0}' ../processed_data/031_find_optimal_PCs/gene_loc.txt > ../processed_data/031_find_optimal_PCs/chr20.gene_loc.txt
subl 031_find_optimal_PCs.R
Rscript 031_find_optimal_PCs.R


# 2016/05/10
# added the following papers to hcasmc binder: 
# 1. uk10k
# 2. dgn
# 3. geuvadis
# 4. immune cell eqtl, fairfax (2012)
# 5. beta cell eqtl, Nica (2013)
# 6. peer

# 11:06
# emailed Milos for an update on the status of RNAseq remap

# 11:07
# started writing the whole-genome sequencing analysis section of the supplement


# 2016/05/11 (Wed)
# 0:18am
# I plan to re-run joint genotype calling and variant recalibration.
# in the variant recalibration step, InbreedingCoeff should be used only when 
# no closely related samples are present. 
# need to check this using PLINK or BEAGLE. 


# calculate IBD using PLINK:
mkdir ../processed_data/160511_calc_ibd_plink
subl 160511_calc_ibd_plink.sh
mkdir ../figures/160511_analyze_ibd
subl 160511_analyze_ibd.R



# 11:14am
# Lab meeting
# Stephen suggested to filter on MAF and genotype quality. 
# And also plot the MAF distribution of pruned set of variants.
# Logic is that if the variants are rare then people are more likely to 
# have IBD2. 


# 11:57am 
# plotting the allele frequency distribution of pruned marker set: 
mkdir ../processed_data/160511_calc_freq
subl 160511_calc_freq.sh
mkdir ../figures/160511_calc_freq/
subl 160511_plot_allele_frequency_of_pruned_marker_set.R 

# Figure ../figures/160511_calc_freq/maf_dist.pdf shows that majority of 
# variants have maf < 0.05. This increases sharing of IBS2 (think the extreme 
# case of maf=0.0). I thus filtered on maf=0.05 and geno=0.1 and recalculated 
# IBD probability
subl 160511_calc_ibd_plink.sh

# I plotted pi_hat and p(IBD=1)
subl 160511_analyze_ibd.R
# Figure ../figures/160511_analyze_ibd/pi_hat.pdf shows that besides three contaminated samples
# 1848, 1868 and 24635, other sample pairs have pi_hat less than 0.125 (more remote
# than 3rd degree relatives). Points very high are duplicates as labeled.

# To understand why contaminated samples have high IBD sharing:
# Figure ../figures/160511_analyze_ibd/Z1.pdf shows that contaminated samples have increases IBD1
# sharing due to excessive heterozygosity.

# It was not clear whether sample 24156 are contaminated. Figure ../figures/160511_analyze_ibd/24156.pdf
# shows that 24156 is unrelated to majority of samples. It is unlikely to have been contaminated.


# 2:19pm
# With IBD information, it is worth recalling and recalibrating genotypes 
# using individuals besides contaminated samples and duplicates
# contaminated samples:
# 1848
# 1858
# 24635

# duplicate pairs: 
# 1346 (180G) and CA1346 (177G)
# 150328 (173G) and 59386145 (173G)
# 2102 (220G) and 2105 (155G)
# 2109 (182G) and CA1508 (175G)
# 289727 (179G) and 2999 (168G)
# 313605 (174G) and 317155 (167G)

# I remove ones with smaller number of mapped reads as indicated 
# in size of the bam files.
# duplicate to remove:
open ../figures/160511_analyze_ibd/pi_hat.pdf
# CA1346
# 150328
# 2105 
# CA1508
# 2999
# 317155 


# 3:38pm
# re-run joint genotyping and recalibration: 
subl genotype_gvcfs.py
mkdir ../data/joint2
cp /mnt/data/WGS_HCASMC/sample_list.txt ../data/joint2/sample_list.txt
subl ../data/joint2/sample_list.txt
python genotype_gvcfs.py ../data/joint2 sample_list.txt 2> log/genotype_gvcfs.20160511.log 
# killed at chromosome 19...
# see 2016/05/12 2:52pm

# 5:21pm 
# added the matrix eQTL paper to hcasmc binder
# I need to understand how to orthogolize out covariates 
# Seber and Lee pg54 is a good read 


# 5:53pm
# BEAGLE genotype refinement sample code: 
wget https://faculty.washington.edu/browning/beagle/run.beagle.03May16.862.example 
mv run.beagle.03May16.862.example beagle_example.sh
subl beagle_example.sh

# 10:50pm
# testing the beagle's conform-gt module:
# I flipped the alleles of rs138720731, and conform-gt successfully flipped back
# I also flipped the strand of rs73387790, and conform-gt reports OPPOSITE_STRAND
# For alleles not in the reference panel, conform-gt reports REMOVED and NOT_IN_REFERENCE
# Note that conform-gt reports FAIL for some correct SNPs likely because 
# the evidence for correct strand is inconclusive. 
# Since conform-gt incorrectly removes these "FAILED" variants, I need to construct vcf manually


# 11:16pm 
# download the beagle reference file: 
mkdir /srv/persistent/bliu2/tools/beagle/reference
cd /srv/persistent/bliu2/tools/beagle/reference
for i in $(seq 1 22) X; do
echo $i
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr$i.1kg.phase3.v5a.bref
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr$i.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr$i.1kg.phase3.v5a.vcf.gz.tbi
done 

# download the genetic map: 
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

# 11:21pm
# read the BEAGLE genotype imputation paper (AJHJ 2016)
# takeaway is beagle, impute2 and minimac3 have similar accuracy but 
# beagle is faster


# 2016/05/12
# 9:38am
# remap dASE samples using Milo's pipeline, making sure STAR versions are the same. 
# download STAR 2.4.0j:
cd /srv/persistent/bliu2/tools/
wget https://github.com/alexdobin/STAR/archive/STAR_2.4.0j.tar.gz
tar -xzf STAR_2.4.0j.tar.gz

# Milos will help map the dASE samples on valkyr: 
mkdir /home/diskstation/RNAseq/dase
screen -S transfer_to_valk
rsync -vzh /srv/persistent/bliu2/dase/data/RNAseq_HCASMC/fastq/*.fastq.gz bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/dase

# Merge fastq files for CA2305:
screen -S merge_fastq
cd /srv/persistent/bliu2/dase/data/RNAseq_CA2305_NextSeq_20150512/fastq
for file in FBS2_S4 FBS4_S5 FBS6_S6 SF4_S1 SF5_S2 SF6_S3; do 
	for read in R1 R2; do 
	echo sample: $file
	echo reads: $read 
	zcat ${file}_L00{1,2,3,4}_${read}_001.fastq.gz | gzip > ${file}_merged_${read}_001.fastq.gz &
	done 
done 
# Transfer to valk:
rsync -vzh /srv/persistent/bliu2/dase/data/RNAseq_CA2305_NextSeq_20150512/fastq/*merged*.fastq.gz bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/dase


# 12:36pm 
# Checked Milo's STAR mapping parameters. He used
# --sjdbOverhang 100 but the mate length is 125bp. 
# This is okay. As long as sjdbOverhang > seedSearchStartLmax (default 50bps), 
# STAR is able to map the seed.


# 2:52pm
# genotype_gvcf.py failed presumably due to memory issues. 
# I will rewrite genotype_gvcf.py to process each chromosome separately
subl genotype_gvcfs.sh
dat=$(date '+%Y%m%d')
for chr in $(seq 1 22) X Y M; do 
bash genotype_gvcfs.sh chr$chr 2> log/genotype_gvcfs.$chr.$dat.log &
done

# 5:30pm 
# read STAR aligner paper by Alex Dobin, added to hcasmc binder



# 2016/05/13
# 9:37am 
# concatenate raw_variants.chr*.vcf.gz
screen -r concat
bcftools concat -Ov -o raw_variants.vcf raw_variants.chr{M,{1..22},X,Y}.vcf
mkdir by_chrom
mv raw_variants.chr*.vcf* by_chrom


# (pilot) recalibrate SNPs on chr22: 
screen -S recalibrate_variants
subl recalibrate_SNP.sh
bash recalibrate_SNP.sh raw_variants.chr22.vcf


# recalibrate SNPs on whole genome:
# I added tranches at 0.1 increments in 99.9-99.0. It seems like 99.5 is the sweet spot
bash recalibrate_SNP.sh raw_variants.vcf


# recalibrate INDELs: 
screen -S recalibrate_INDEL
subl recalibrate_INDEL.sh 
bash recalibrate_INDEL.sh recalibrated_snps_raw_indels.vcf 2> log/recalibrate_INDEL.160513.sh

# plot sensitivity vs minVOSLod score: 
mkdir ../figures/160513_plot_sensitivity/
subl plot_sensitivity.R
Rscript plot_sensitivity.R
# conclusion: seems like 98% sensitivity (elbow) is a good cutoff. 

# run ApplyRecalibration (commented out VariantRecalibrator code):
bash recalibrate_INDEL.sh recalibrated_snps_raw_indels.vcf 2> log/recalibrate_INDEL.160514.sh

# remove intermediate files:
wd=../data/joint2/
rm $wd/recalibrated_snps_raw_indels.vcf
rm $wd/recalibrated_snps_raw_indels.vcf.idx
mkdir $wd/variant_recalibration
mv $wd/recalibrate_SNP* $wd/variant_recalibration
mv $wd/recalibrate_INDEL* $wd/variant_recalibration


# 4:57pm
# checking possible sample contamination
mkdir ../processed_data/160513_contamination
verifyBamID --vcf /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --bam /mnt/data/WGS_HCASMC/1020301/recal_reads.bam --chip-none --precise --verbose --minAF 0 --minCallRate 0 --out ../processed_data/160513_contamination/chr22
zcat ../processed_data/dna_contamination/chr20.AF.EUR.2.vcf.gz | head -n10000 > ../processed_data/160513_contamination/chr20.AF.EUR.2.head10000.vcf
verifyBamID --vcf ../processed_data/160513_contamination/chr20.AF.EUR.2.head10000.vcf --bam /mnt/data/WGS_HCASMC/1020301/recal_reads.bam --chip-none --precise --verbose --minAF 0 --minCallRate 0 --out ../processed_data/160513_contamination/chr22
# subset bam to chr20
samtools view -h /mnt/data/WGS_HCASMC/1020301/recal_reads.bam chr20:1-100000 | sed 's/chr20/20/' | samtools view -h -b -o ../processed_data/160513_contamination/recal_reads.chr20.bam -
samtools index ../processed_data/160513_contamination/recal_reads.chr20.bam
zcat /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | head -n10000 > ../processed_data/160513_contamination/chr20.vcf
verifyBamID --vcf ../processed_data/160513_contamination/chr20.vcf --bam ../processed_data/160513_contamination/recal_reads.chr20.bam --chip-none --precise --verbose --minAF 0 --minCallRate 0 --out ../processed_data/160513_contamination/chr22
# I decided not to finish this analysis because of the following difficulties
# 1) for admixed samples, MAF cannot be obtained from 1000 genome
# 2) called genotypes may be incorrect


# 6:26pm 
# What is considered SNP by bcftools? 
for type in snps indels mnps other;do
bcftools view -Ov -o /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/by_chrom/raw_variants.chr20.$type.vcf --types $type /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/by_chrom/raw_variants.chr20.vcf &
done
# other means the following for example: 
# chr20   83252   rs6137896       G       C,<*:DEL>


# 2016/05/14:
# 4:46pm
# filter out variant that failed the filters:
screen -S filter_variants
wd=../data/joint2/
bcftools view -Ov -o $wd/recalibrated_variants.pass.vcf -f PASS $wd/recalibrated_variants.vcf


# calculate Ti/Tv ratio:
mkdir ../processed_data/160514_Ts_Tv_ratio
bcftools stats $wd/recalibrated_variants.pass.vcf > ../processed_data/160514_Ts_Tv_ratio/stats.txt
plot-vcfstats --no-PDF -p ../processed_data/160514_Ts_Tv_ratio/stats ../processed_data/160514_Ts_Tv_ratio/stats.txt


# what proportion of variants are biallelic SNPs: 
screen -S count_variant_by_type
mkdir ../processed_data/160514_variant_count_by_type
mkdir ../figures/160514_variant_count_by_type
cat ../processed_data/160514_Ts_Tv_ratio/stats.txt | grep "^SN" | awk 'BEGIN {FS="\t"} {print $0}' > ../processed_data/160514_variant_count_by_type/count.txt
subl 160514_plot_variant_count_by_type.R


# filter for biallelic snps:
screen -S filter_biallelic_SNP
bcftools view -m2 -M2 -v snps -Ov -o $wd/recalibrated_biallelic_SNP.pass.vcf $wd/recalibrated_variants.pass.vcf 


# 7:00pm
# impute with BEAGLE without reference panels:
screen -S beagle_no_ref
wd=../data/joint2/
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java
beagle=/srv/persistent/bliu2/tools/beagle/beagle.03May16.862.jar
for i in $(seq 1 22);do
$java -Xmx4096m -jar $beagle nthreads=2 chrom=chr$i gl=$wd/recalibrated_biallelic_SNP.pass.vcf out=$wd/recalibrated_biallelic_SNP.beagle.chr$i &
done
# index each vcf:
for i in $(seq 1 22); do
tabix -p vcf recalibrated_biallelic_SNP.beagle.chr$i.vcf.gz
done 
# concatenate all chromosomes:
bcftools concat -Ov -o recalibrated_biallelic_SNP.beagle.vcf recalibrated_biallelic_SNP.beagle.chr{1..22}.vcf.gz
# I made sure the beagle output is complete by comparing the CHROM and POS fields of 
# recalibrated_biallelic_SNP.beagle.vcf and recalibrated_biallelic_SNP.pass.vcf
# move beagle intermediate files to folder: 
mkdir beagle_no_ref
mv recalibrated_biallelic_SNP.beagle.chr* beagle_no_ref


# 2016/05/15:
# change chromosome names from "chr*" to "*":
wd=../data/joint2/
for i in $(seq 1 22) X Y M; do 
echo "chr$i $i" >> $wd/old_to_new_chrom_name.txt
done
screen -S change_chrom_name
bcftools annotate --rename-chrs $wd/old_to_new_chrom_name.txt -Ov -o $wd/recalibrated_biallelic_SNP.pass.GRCh37.vcf $wd/recalibrated_biallelic_SNP.pass.vcf
# beagle with reference panel:
screen -S beagle_with_ref
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java
beagle=/srv/persistent/bliu2/tools/beagle/beagle.03May16.862.jar
reference_dir=/srv/persistent/bliu2/tools/beagle/reference
for i in $(seq 1 22);do
i=22 # test on chr22
$java -Xmx32g -jar $beagle nthreads=24 chrom=$i ref=$reference_dir/chr$i.1kg.phase3.v5a.vcf.gz map=$reference_dir/plink.chr$i.GRCh37.map impute=false gl=$wd/recalibrated_biallelic_SNP.pass.GRCh37.vcf out=$wd/recalibrated_biallelic_SNP.beagle_1kg.chr$i > $wd/recalibrated_biallelic_SNP.beagle_1kg.chr$i.log2 &
done


# 11:32am 
# beagle output QC:
# to plot the dosage R2 as a function of chromosome, allele frequency
# also to make histogram of dosage R2
screen -S beagle_QC
wd=../data/joint2/
mkdir ../processed_data/160515_beagle_QC
mkdir ../figures/160515_beagle_QC/
cat $wd/recalibrated_biallelic_SNP.beagle.vcf | grep -v "^#" | awk 'BEGIN {FS="\t|;|="; OFS="\t"; print "CHROM","POS","AR2","DR2","AF"} {print $1,$2,$9,$11,$13}' >  ../processed_data/160515_beagle_QC/recalibrated_biallelic_SNP.r2.tsv
subl beagle_QC.R
Rscript beagle_QC.R
# reference on how to calculate genotype R-squared: http://ingenoveritas.net/compare-true-and-imputed-genotypes-by-calculating-r-squared/


# 5:45pm 
# update sample names:
wd=../data/joint2/
cd $wd
echo "CA1401 1401" >> old_to_new_sample_name.txt
echo "2102 2105" >> old_to_new_sample_name.txt
echo "2109 1508" >> old_to_new_sample_name.txt
echo "289727 2999" >> old_to_new_sample_name.txt
echo "313605 317155" >> old_to_new_sample_name.txt
screen -S rename
bcftools reheader -s old_to_new_sample_name.txt -o recalibrated_biallelic_SNP.beagle.rename.vcf recalibrated_biallelic_SNP.beagle.vcf

# filter for variants with dosage R2 >= 0.8:
bcftools view -e 'INFO/DR2<0.8' -o recalibrated_biallelic_SNP.beagle.rename.dr2.vcf recalibrated_biallelic_SNP.beagle.rename.vcf
 
# subset for Caucasian individuals to apply HWE filter: 
subl caucasian_individual_for_hwe.R


# make sample list
# list only contains unique sample names sorted alphanumerically for the 52 samples in the working set
subl make_sample_list.R

# extract dosage field:
# the output column order will be the same as that in sample_list.txt
wd=../data/joint2/
cd /srv/persistent/bliu2/HCASMC_eQTL/scripts
mkdir ../processed_data/160515_dosage
screen -S extract_DS
bcftools query -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%DS]\n' -o ../processed_data/160515_dosage/dosage.tsv $wd/recalibrated_biallelic_SNP.beagle.rename.vcf &
bcftools query -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' -o ../processed_data/160515_dosage/genotype.tsv $wd/recalibrated_biallelic_SNP.beagle.rename.vcf &
bcftools query -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%GP]\n' -o ../processed_data/160515_dosage/genotype_probability.tsv $wd/recalibrated_biallelic_SNP.beagle.rename.vcf &




# 7:49pm
# archive some files:
screen -S archive
bgzip raw_variants.vcf &
bgzip recalibrated_variants.vcf &
bgzip recalibrated_variants.pass.vcf &
bgzip recalibrated_biallelic_SNP.pass.vcf &
bgzip recalibrated_biallelic_SNP.beagle.vcf &



# 2016/05/16
# Milos used STAR v2.5.1 instead of v2.4.0. So I need to remap. 
# The three samples I need are FBS2_S4_merged_R1_001.fastq.gz (most reads among replicates)
# S7_run0002_lane5_index7_1.fastq.gz (20805), pS17_1.fastq.gz (9052004)
# on valk: 
cd /home/diskstation/RNAseq/dase
cp commands commands2
vim commands2 
# changed STAR to $STAR for pass2
screen -S STAR
bash commands2
# I don't have permission to write...


# 12:00pm
# detect RNAseq experssion outliers: 
mkdir ../figures/160516_detect_expression_outlier/
subl 160516_detect_expression_outlier.R
# 2135, 2305 and 9070202 are outliers
# 9070202 have low mapping rate (23%) so should use the remapped reads
# 2135 has double peaked bioanalyzer result
# 2305 is sequenced on a NextSeq separately from all other samples.


# does omitting covariate decrease power?
mkdir ../figures/160516_sim_study_on_covariates/
subl 160516_sim_study_on_covariates.R
# conclusion: omitting covariates will decrease power.


# 5:48pm
# prepare matrix eQTL genotype: 
wd=../data/joint2
dir1=../processed_data/160516_genotype
mkdir $dir1
bcftools query -H -e 'INFO/DR2<0.8' -t chr22 -S $wd/sample_list.txt -f '%CHROM\_%POS\_%REF\_%ALT[\t%DS]\n' -o $dir1/chr22.gneotype.tmp $wd/recalibrated_biallelic_SNP.beagle.rename.vcf
sed -e "s/# \[1\]CHROM_\[2\]POS_\[3\]REF_\[4\]ALT/id/" -e "s/\[[[:digit:]]\+\]//g" -e "s/:DS//g" -e "s/chr//" $dir1/chr22.gneotype.tmp > $dir1/chr22.genotype.txt

# filter for genotype with maf>=0.05:
subl 160516_subset_genotype_by_maf.R 
Rscript 160516_subset_genotype_by_maf.R $dir1/chr22.genotype.txt $dir1/chr22.genotype.maf.txt

# prepare snp location:
subl 160516_gen_snps_loc.R
Rscript 160516_gen_snps_loc.R $dir1/chr22.genotype.maf.txt $dir1/chr22.genotype_loc.maf.txt

# prepare gene expression:
# gene expression is already prepared 


# 7:00pm
# on 16/05/16 12:00pm we determined that there are 3 outliers, 
# here we determine whether including them will decrease power.
# prepare genotype and expression files without the 3 outliers, which 
# are columns 30, 34, and 50
cut -f-29,31-33,35-49,51- $dir1/chr22.genotype.maf.txt > $dir1/chr22.genotype.maf.cut.txt
cut -f-29,31-33,35-49,51- ../processed_data/031_prepare_matrix_eQTL_expression/expression.txt > $dir1/expression.cut.txt
screen -S matrixeqtl
mkdir ../figures/160516_matrix_eQTL/
subl 160516_matrix_eQTL.R $dir1/chr22.genotype.maf.txt $dir1/chr22.genotype_loc.maf.txt ../processed_data/031_prepare_matrix_eQTL_expression/expression.txt ../processed_data/031_gen_gene_loc/gene_loc.txt "" $dir1


# 2016/05/17
# 6:41pm
# get a clean set of gene expression data:
mkdir ../data/rnaseq2
# saved ../processed_data/rna_wgs_match.reduced_050616.xlsx into txt file 
# use vim to turn all ^M into \r
subl 160517_get_rnaseq_data.sh

# sort and index some bam and sam files: 
samtools sort -o 2305/Aligned.out.sorted.bam -O bam -@8 2305/Aligned.out.sam &
samtools sort -o 9070202/Aligned.out.sorted.bam -O bam -@8 9070202/Aligned.out.bam &
samtools sort -o 9052004/Aligned.out.sorted.bam -O bam -@8 9052004/Aligned.out.bam &
samtools sort -o 20805/Aligned.out.sorted.bam -O bam -@8 20805/Aligned.out.sam &
samtools index 2305/Aligned.out.sorted.bam & 
samtools index 9070202/Aligned.out.sorted.bam & 
samtools index 9052004/Aligned.out.sorted.bam & 
samtools index 20805/Aligned.out.sorted.bam & 
rm 2305/Aligned.out.sam 9070202/Aligned.out.bam 9052004/Aligned.out.bam 20805/Aligned.out.sam

# cehck that all files are intact:
samtools quickcheck */Aligned.out.sorted.bam
# all files are intact. 


# 11:35pm
# create sequence dictionary:
CreateSequenceDictionary=/software/picard-tools/1.92/CreateSequenceDictionary.jar
java -jar $CreateSequenceDictionary R=/srv/persistent/bliu2/shared/genomes/hg19/hg19.fa O=/srv/persistent/bliu2/shared/genomes/hg19/hg19.dict

# add read group:
screen -S readgroup
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
AddOrReplaceReadGroups=/software/picard-tools/1.92/AddOrReplaceReadGroups.jar
samples=($(ls))
n=0
for sample in ${samples[@]}; do 
echo $sample
n=$((n+1))
if [[ n -gt 30 ]];then
	wait
	n=0
fi 
java -Xmx2g -jar $AddOrReplaceReadGroups I=$sample/Aligned.out.sorted.bam O=$sample/Aligned.out.sorted.rg.bam RGID=$sample RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample 2> AddOrReplaceReadGroups.$sample.log &
done
wait

n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 30 ]];then
	wait
	n=0
fi 
echo $sample
samtools index $sample/Aligned.out.sorted.rg.bam &
done


# filter for uniquely mapped reads: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
samples=($(ls -d */))
screen -S uniquely_mapped
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 15 ]];then
	wait
	n=0
fi
echo $sample
samtools view -h -q 255 -b -o $sample/Aligned.out.sorted.rg.uniq.bam $sample/Aligned.out.sorted.rg.bam &
done

n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 20 ]];then
	wait
	n=0
fi 
echo $sample
samtools index $sample/Aligned.out.sorted.rg.uniq.bam &
done

# mark duplicates: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
screen -S mark_dup
samples=($(ls -d */))
MarkDuplicates=/software/picard-tools/1.92/MarkDuplicates.jar 
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 15 ]];then
	wait
	n=0
fi
echo $sample
java -Xmx2g -jar $MarkDuplicates I=$sample/Aligned.out.sorted.rg.uniq.bam O=$sample/Aligned.out.sorted.rg.uniq.dup.bam  M=$sample/marked_dup_metrics.txt 2> MarkDuplicates.${sample///}.log &
done

n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 30 ]];then
	wait
	n=0
fi 
echo $sample
samtools index $sample/Aligned.out.sorted.rg.uniq.dup.bam &
done

# make sample file for RNA-seQC:
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
ls -d */ > sample_file.tmp
cat sample_file.tmp | sed "s:/::" | awk 'BEGIN {OFS="\t"; print "Sample ID","Bam File","Notes"} {print $1,$1"/Aligned.out.sorted.rg.uniq.dup.bam",$1}' > sample_file.txt
rm sample_file.tmp


# get gtex v6p gencode gene models: 
# this annotation collapses exons into genes:
cd /srv/persistent/bliu2/shared/annotation/
mkdir gtex
ln -s /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/reference_files/gencode.v19.genes.v6p.patched_contigs.gtf.gz . 
gunzip -c gencode.v19.genes.v6p.patched_contigs.gtf.gz > gencode.v19.genes.v6p.patched_contigs.gtf
cat gencode.v19.genes.v6p.patched_contigs.gtf | sed -e "/^[^#]/ s/^/chr/" -e "s/MT/M/" > gencode.v19.genes.v6p.hg19.gtf


# run RNA-seQC: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments
screen -S rnaseqc
rnaseqc=/srv/persistent/bliu2/tools/RNA-SeQC_v1.1.8.jar
gencode19=/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf
hg19=/srv/persistent/bliu2/shared/genomes/hg19/hg19.fa
rRNA=/srv/persistent/bliu2/shared/genomes/rRNA/human_all_rRNA.fasta
# java -jar $rnaseqc -n 1000 -s sample_file.txt -t $gencode19 -r $hg19 -o report -noDoC -strictMode

samples=($(ls -d */))
n=0
for sample in ${samples[@]}; do 
n=$((n+1))
if [[ n -gt 20 ]];then
	wait
	n=0
fi
sample=${sample///} # remove the backslash

if [[ $(grep "Finished Successfully" rnaseqc.$sample.log) == "" ]]; then 
echo $sample
java -jar $rnaseqc -n 1000 -s "$sample|$sample/Aligned.out.sorted.rg.uniq.dup.bam|$sample" -t $gencode19 -r $hg19 -o $sample/report -noDoC -strictMode > rnaseqc.$sample.log &
fi 
done

# 16/05/19:
# combine all RPKMs:
mkdir ../processed_data/160519_rpkm/
wd=/srv/persistent/bliu2/HCASMC_eQTL/data/rnaseq2/alignments

tail -n +3 $wd/$sample/report/genes.rpkm.gct | cut -f1-2 > ../processed_data/160519_rpkm/combined.rpkm
samples=($(ls -d $wd/*/))

for sample in ${samples[@]};do
sample=$(basename $sample)
sample=${sample///}
echo $sample
tail -n +3 $wd/$sample/report/genes.rpkm.gct | cut -f3 > ../processed_data/160519_rpkm/$sample.rpkm.tmp
cp ../processed_data/160519_rpkm/combined.rpkm ../processed_data/160519_rpkm/combined.rpkm.tmp
paste -d "\t" ../processed_data/160519_rpkm/combined.rpkm.tmp ../processed_data/160519_rpkm/$sample.rpkm.tmp > ../processed_data/160519_rpkm/combined.rpkm
done
rm ../processed_data/160519_rpkm/*.tmp


# calculate the sample-sample correlation:
mkdir ../figures/160519_rpkm
subl 160519_calc_sample_correlation.R


# 5:01pm 
# convert VCF to plink BED file:
mkdir ../processed_data/160519_genotype_PCA
plink --vcf /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.vcf --keep-allele-order --make-bed --out ../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename
plink --vcf /srv/persistent/bliu2/HCASMC_eQTL/data/joint2/recalibrated_biallelic_SNP.beagle.rename.dr2.vcf --keep-allele-order --make-bed --out ../processed_data/160519_genotype_PCA/recalibrated_biallelic_SNP.beagle.rename.dr2

# find genotype PCs: 
# Bruna shared script to call genotype PCs through slack:
mkdir ../figures/160519_genotype_PCA
subl 160519_genotype_PCA.R
# seems that 1020301 is an outlier? It has low heterozygosity rate (from the plink/seq analysis)

# 16/05/20
# 11:52am
# 


# look at WGS unique mapping rate for 1020301. 


# 3:07pm
# On 16/05/19 I showed that sample 9052004 (Stanford 2nd round sequencing) is an outlier. 
# Here we analyze the dASE version of 9052004. 
mv ../data/rnaseq2/alignments/9052004 ../data/rnaseq2/alignments/.9052004 # hide the bad sample
dst=../data/rnaseq2/alignments/
rsync -azvh bosh@valkyr.stanford.edu:/home/diskstation/RNAseq/dase/pS17_1.fastq.gz_pS17_2.fastq.gz/Pass2/{Aligned.out.sam,Log.final.out,Log.out,Log.progress.out} $dst/9052004 &
subl rerun_rnaseqc_for_9052004_dase.sh
screen -S rerun_for_9042004
bash rerun_rnaseqc_for_9052004_dase.sh

wd=../data/rnaseq2/alignments/
mkdir ../processed_data/160520_rpkm/
tail -n +3 $wd/1020301/report/genes.rpkm.gct | cut -f1-2 > ../processed_data/160520_rpkm/combined.rpkm
samples=($(ls -d $wd/*/))
for sample in ${samples[@]};do
sample=$(basename $sample)
sample=${sample///}
echo $sample
tail -n +3 $wd/$sample/report/genes.rpkm.gct | cut -f3 > ../processed_data/160520_rpkm/$sample.rpkm.tmp
cp ../processed_data/160520_rpkm/combined.rpkm ../processed_data/160520_rpkm/combined.rpkm.tmp
paste -d "\t" ../processed_data/160520_rpkm/combined.rpkm.tmp ../processed_data/160520_rpkm/$sample.rpkm.tmp > ../processed_data/160520_rpkm/combined.rpkm
done
rm ../processed_data/160520_rpkm/*.tmp


# calculate the sample-sample correlation:
mkdir ../figures/160520_rpkm
subl 160520_calc_sample_correlation.R
Rscript 160520_calc_sample_correlation.R
