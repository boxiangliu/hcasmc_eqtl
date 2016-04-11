from utils import *

GATK="/usr/bin/GenomeAnalysisTK.jar"
HG19 = '/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
GRCH37 = '/srv/persistent/bliu2/shared/genomes/GRCh37/hs37d5.fa'

def calcGenotypePosteriors(in_vcf, out_vcf, supporting="/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.concat.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz", dry_run = False):
	cmd = "java -Xmx4g -jar %s -T CalculateGenotypePosteriors -R %s --supporting %s -V %s -o %s"%(GATK, GRCH37, supporting, in_vcf, out_vcf)
	call(cmd, dry_run)
	print "[ calcGenotypePosteriors ] finished"

def variantFiltrationByGQ(in_vcf, out_vcf, GQ = 20.0, dry_run = False):
	cmd = 'java -Xmx4g -jar %s -T VariantFiltration -R %s -V %s -G_filter "GQ < %s" -G_filterName lowGQ -o %s'%(GATK, GRCH37, in_vcf, GQ, out_vcf)
	call(cmd, dry_run)
	print "[ variantFiltrationByGQ ] finished"

if __name__ == '__main__':
	dry_run = False
	
	printDateTime()
	wd = sys.argv[1]
	in_vcf = sys.argv[2]
	setwd(wd)

	CGP_vcf = in_vcf.replace('vcf', 'postCGP.vcf').replace('.gz', '')
	# calcGenotypePosteriors(in_vcf, CGP_vcf, dry_run= dry_run)

	filtered_vcf = CGP_vcf.replace('vcf', 'GQfiltered.vcf')
	# variantFiltrationByGQ(CGP_vcf, filtered_vcf, dry_run = dry_run)

	compressVcf(sample_dir = './', sample = filtered_vcf)
	indexVcf(sample_dir = './', sample = filtered_vcf + '.gz')
	reportDone()