#!/usr/bin/env python 
import gzip,sys
sys.path.append('/srv/persistent/bliu2/HCASMC_eQTL/scripts')
from utils import *
debug = False
# read VCF file:

def readVCF(vcf):
	snps = []
	with gzip.open(vcf, 'r') as f: 
		for line in f: 
			line = line.strip()
			if line.startswith("##"): continue
			split_line = line.split('\t')
			if line.startswith("#"):
				chrom_idx = split_line.index('#CHROM')
				pos_idx = split_line.index('POS')
				ref_idx = split_line.index('REF')
				alt_idx = split_line.index('ALT')
			else:
				chrom = split_line[chrom_idx]
				pos = split_line[pos_idx]
				ref = split_line[ref_idx]
				alt = split_line[alt_idx]
				snps.append([chrom, pos, ref, alt])
	return snps


# write to SNP file: 
def writeSNPfile(snps):
	last_chrom = ""
	for snp in snps: 
		chrom = snp[0]
		if chrom != last_chrom: # open a new file for each chromosome. 
			try: 
				out.close()
			except NameError:
				pass  
			out = gzip.open("%s.snps.txt.gz"%chrom, 'w')
			print "writing chr%s snps"%chrom
		pos = snp[1]
		ref = snp[2]
		alt = snp[3]
		out.write(pos + '\t' + ref + '\t' + alt + "\n")

		last_chrom = chrom
	
	out.close()

vcf=sys.argv[1]
outdir=sys.argv[2]
setwd(outdir)
print '[ input VCF ]', vcf
snps = readVCF(vcf)
writeSNPfile(snps)
reportDone()



