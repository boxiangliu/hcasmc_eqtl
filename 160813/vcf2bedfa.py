# Convert VCF file to bed-fasta pair: 
# 
# Author: Boxiang Liu
# Email: bliu2@stanford.edu
# 08/13/2016
# 
# Input arguments: 
# vcf_file: the vcf file to be converted to bed-fasta pairs
# sample: the sample name. It should be a genotype column of the VCF file
# out_prefix: prefix to the output. 
# 
# Output: 
# two files will be generated. <output_prefix>.bed and <output_prefix>.fa. 
 
import vcf,sys

# command args: 
vcf_file=sys.argv[1]
sample=sys.argv[2]
out_prefix=sys.argv[3]

# example inputs: 
# vcf_file='/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/recalibrated_biallelic_variants.beagle.rename.dr2.hwe.maf.head2000.vcf'
# sample='20805'
# out_prefix='/srv/persistent/bliu2/HCASMC_eQTL/data/joint3/bed_fa/20805'


# format output file names: 
bed_file=out_prefix+'.bed'
fa_file=out_prefix+'.fa'


# get vcf object:
vcf_reader = vcf.Reader(open(vcf_file,'r'))
vcf_reader.samples


# iterate over vcf file:
with open(bed_file,'w') as bed, open(fa_file,'w') as fa:
	n_inline=0
	for record in vcf_reader:
		n_inline=n_inline+1

		# get chrom, start, end
		chrom=record.CHROM
		start=record.POS
		ref_length=len(record.REF)
		end=start+ref_length-1

		# cast all variants into a list: 
		variants=[record.REF]+[str(x) for x in record.ALT]

		# cast genotype into a sorted list: 
		geno=record.genotype(sample)['GT']
		geno=geno.replace('|',"/").split('/')
		geno=[int(x) for x in geno]
		geno.sort()

		# write to bed file: 
		bed_out="\t".join([str(chrom),str(start-1),str(end)])
		bed.write(bed_out+'\n')

		# write to fasta file:
		fa_out=">"+str(chrom)+":"+str(start-1)+"-"+str(end)+'\n'+variants[geno[-1]]
		fa.write(fa_out+'\n')

		# progress report: 
		if n_inline%100000==0:
			sys.stderr.write(str(n_inline)+' lines\n')
