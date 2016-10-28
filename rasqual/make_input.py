####
# python form_input.py <window_size>
#### 
import sys
sys.path.insert(1,'/srv/persistent/bliu2/tools/anaconda/lib/python2.7/site-packages/')
import pysam,vcf,re

# functions: 
def count_fsnp_and_rsnp(chrom,window_start,window_end,exons):
	# count number of test and feature snps: 
	n_test_snp=0
	n_feat_snp=0
	for record in vcf_reader.fetch(chrom,window_start-1,window_end):
		n_test_snp+=1
		is_fsnp=reduce(lambda x,y: x or y, [exon[0]<=record.POS<=exon[1] for exon in exons])
		if is_fsnp:
			n_feat_snp+=1
	return n_test_snp,n_feat_snp

def make_output(gene_id,gene_name,chrom,window_start,window_end,n_test_snp,n_feat_snp,exons):
	# make output: 
	exon_start_positions=','.join([str(exon[0]) for exon in exons])
	exon_end_positions=','.join([str(exon[1]) for exon in exons])
	output="\t".join([gene_id,gene_name,chrom+":"+str(window_start)+"-"+str(window_end),str(n_test_snp),str(n_feat_snp),exon_start_positions,exon_end_positions])+'\n'
	return output


# cmd args: 
vcf_file=sys.argv[1]
window_size=int(sys.argv[2])


# read VCF file: 
vcf_reader = vcf.Reader(open(vcf_file, 'r'))


# regex to find gene name: 
p_gene_id=re.compile('(?<=gene_id\s").+?(?=";)')
p_gene_name=re.compile('(?<=gene_name\s").+?(?=";)')


for line in sys.stdin:
	if "##" in line: # skip headers
		continue
	else: 
		split_line=line.strip().split('\t')

	if split_line[2]=='gene': # gene lines

		#### output: 
		try:
			# report: 
			sys.stderr.write(gene_id+' '+gene_name+' has '+str(len(exons))+' exon(s)\n')

			# count number of test and feature snps: 
			[n_test_snp,n_feat_snp]=count_fsnp_and_rsnp(chrom, window_start, window_end, exons)
			output=make_output(gene_id, gene_name, chrom, window_start, window_end, n_test_snp, n_feat_snp, exons)
			sys.stdout.write(output)

		except NameError: 
			pass

		except ValueError: 
			sys.stderr.write(chrom+' not found in '+vcf_file+'\n')

		#### update: 
		# gene id and name: 
		gene_id=p_gene_id.search(split_line[8]).group()
		gene_name=p_gene_name.search(split_line[8]).group()
		
		# retrieve TSS:
		if split_line[6]=='+': # plus strand
			tss=int(split_line[3])
		elif split_line[6]=='-': # minus strand
			tss=int(split_line[4])
		else: # missing strand information
			tss=int(split_line[3])
			sys.stderr.write('missing strand information; assuming + strand for TSS.')

		# chromosome and testing window start and end:  
		chrom=split_line[0]
		window_start=max(tss-window_size,1)
		window_end=tss+window_size


		# # back up some variables for next loop: 
		# gene_id_bak=gene_id
		# gene_name_bak=gene_name
		# chrom_bak=chrom
		# window_start_bak=window_start
		# window_end_bak=window_end


		# initialize:
		exons=[]


	elif split_line[2]=='exon': # exon lines 
		exon_start=int(split_line[3])
		exon_end=int(split_line[4])
		exons.append((exon_start,exon_end))

	else: pass 

# final output:
try:
	sys.stderr.write(gene_id+' '+gene_name+' has '+str(len(exons))+' exon(s)\n')
	[n_test_snp,n_feat_snp]=count_fsnp_and_rsnp(chrom, window_start, window_end, exons)
	output=make_output(gene_id, gene_name, chrom, window_start, window_end, n_test_snp, n_feat_snp, exons)
	sys.stdout.write(output)
except ValueError:
	sys.stderr.write(chrom+' not found in '+vcf_file+'\n')


