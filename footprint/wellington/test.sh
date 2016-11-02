# perform footprinting: 
wellington_footprints.py \
../data/atacseq/fbs/2305/out/peak/idr/optimal_set/2305_ppr.IDR0.1.filt.12-col.2.bed \
../data/atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam \
../processed_data/footprint/wellington/


# estimate DNaseI cleavage bias: 
dnase_bias_estimator.py CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.WellingtonFootprints.FDR.0.01.bed \
/srv/persistent/bliu2/HCASMC_eQTL/data//atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam \
/srv/persistent/bliu2/shared/shared/genomes/hg19/hg19.fa \
/srv/persistent/bliu2/shared/genomes/hg19.chrom.sizes \
CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.WellingtonFootprints.FDR.0.01.bed.dnase_bias


# plot DNase average profile:
dnase_average_profile.py CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.WellingtonFootprints.FDR.0.01.bed \
/srv/persistent/bliu2/HCASMC_eQTL/data//atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam \
CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.average_profile.pdf


# plot DNase average profile (with bias correction):
dnase_average_profile.py CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.WellingtonFootprints.FDR.0.01.bed \
/srv/persistent/bliu2/HCASMC_eQTL/data//atacseq/fbs/2305/out/align/rep1/CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam \
CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.average_profile.bias_correction.pdf \
-bf CA2305-FBS1_S1_concat_R1_001.trim.PE2SE.bam.2305_ppr.IDR0.1.filt.12-col.2.bed.WellingtonFootprints.FDR.0.01.bed.dnase_bias.pickle


