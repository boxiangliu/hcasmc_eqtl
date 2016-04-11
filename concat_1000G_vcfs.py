from utils import *

VCFTOOLS = "/usr/local/bin/vcf-concat"

def concatVcf(vcfs, concatenated_vcf, dry_run = False): 
	cmd = VCFTOOLS
	for vcf in vcfs:
		cmd += " " + vcf
	cmd += " > %s"%concatenated_vcf
	call(cmd, dry_run)
	print "[ concatVcf ] finished."

if __name__ == '__main__':
	dry_run = False
	
	wd = sys.argv[1]
	sample_list = sys.argv[2]
	concatenated_vcf = sys.argv[3]
	setwd(wd)
	
	vcfs = readSampleList(sample_list)

	concatVcf(vcfs, concatenated_vcf, dry_run)
	compressVcf(sample_dir = './', sample = concatenated_vcf)
	indexVcf(sample_dir = './', sample = concatenated_vcf + '.gz')
	reportDone()