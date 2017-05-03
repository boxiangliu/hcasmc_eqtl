from utils import *

GATK="/usr/bin/GenomeAnalysisTK.jar"
HG19 = '/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
GATK_BUNDLE = '/srv/persistent/bliu2/shared/gatk_bundle_2.8_hg19'

def genotypeGVCFs(sample_list, out_vcf = 'raw_variants.vcf', dry_run = False):
	cmd = "java -Xmx24g -jar %s -nt 6 -T GenotypeGVCFs -R %s --dbsnp %s/dbsnp_138.hg19.vcf"%(GATK,HG19,GATK_BUNDLE)
	for sample in sample_list:
		cmd += " --variant %s"%sample
	cmd += ' -o %s'%out_vcf
	call(cmd,dry_run)
	print "[ genotypeGVCFs ] finished."

def calibrateSNPs(vcf, recal_vcf, tranche = 99, dry_run = False):
	cmd = "java -Xmx32g -jar %s -nt 8 -T VariantRecalibrator -R %s -input %s -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 %s/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 %s/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R"%(GATK,HG19,vcf,GATK_BUNDLE,GATK_BUNDLE,GATK_BUNDLE,GATK_BUNDLE)
	call(cmd,dry_run)
	cmd = "java -Xmx32g -jar %s -nt 8 -T ApplyRecalibration -R %s -input %s -mode SNP --ts_filter_level %i -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o %s"%(GATK, HG19, vcf, tranche, recal_vcf)
	call(cmd,dry_run)
	print "[ calibrateSNPs ] finished."

def calibrateINDELs(vcf, recal_vcf, tranche = 99, dry_run = False):
	cmd = 'java -Xmx32g -jar %s -nt 8 -T VariantRecalibrator -R %s -input %s -resource:mills,known=true,training=true,truth=true,prior=12.0 %s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 97.0 -tranche 90.0 --maxGaussians 4 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R'%(GATK, HG19, vcf, GATK_BUNDLE)
	call(cmd,dry_run)
	cmd = 'java -Xmx32g -jar %s -nt 8 -T ApplyRecalibration -R %s -input %s -mode INDEL --ts_filter_level %i -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -o %s'%(GATK, HG19, vcf, tranche, recal_vcf)
	call(cmd,dry_run)
	print "[ calibrateINDELs ] finished."

if __name__ == '__main__':
	dry_run = False
	
	wd = sys.argv[1]
	sample_list_file = sys.argv[2]
	
	setwd(wd)
	
	sample_list = readSampleList(sample_list_file)
	sample_list = ['/mnt/data/WGS_HCASMC/'+sample+'/raw_variants.g.vcf.gz' for sample in sample_list]
	for sample in sample_list: print sample

	joint_vcf = 'raw_variants.vcf'
	genotypeGVCFs(sample_list, joint_vcf, dry_run)
	
	recal_snp_raw_indel_vcf = 'recalibrated_snps_raw_indels.vcf'
	calibrateSNPs(joint_vcf, recal_snp_raw_indel_vcf, 99, dry_run)
	
	recal_vcf = 'recalibrated_variants.vcf'
	calibrateINDELs(recal_snp_raw_indel_vcf, recal_vcf, 99, dry_run)

	compressVcf(sample_dir = './', sample = recal_vcf)
	indexVcf(sample_dir = './', sample = recal_vcf + '.gz')
	reportDone()