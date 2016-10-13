####
# python calc_gcc.py genome.fa annotation.gtf <gene|exon>
#### 
from pyfaidx import Fasta
import sys
from Bio.SeqUtils import GC
genome_file=sys.argv[1]
annotation_file=sys.argv[2]
genome_file='/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa'
annotation_file='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf'

genome = Fasta(genome_file)

if sys.argv[3]=='gene':
	with open(annotation_file,'r') as f:
		for line in f:
			if "##" in line: continue
			split_line=line.strip().split('\t')
			if split_line[2]=='gene':
				chrom=split_line[0]
				start=int(split_line[3])-1
				end=int(split_line[4])
				gene_id=split_line[8].split(';')[0].replace('gene_id ','').replace('"','')
				gcc=GC(genome[chrom][start:end].seq)
				sys.stdout.write("%s\t%s\n"%(gene_id,gcc))

elif sys.argv[3]=='exon':
	with open(annotation_file,'r') as f:
		for line in f:
			if "##" in line: continue

			split_line=line.strip().split('\t')
			if split_line[2]=='gene':
				gene_id=split_line[8].split(';')[0].replace('gene_id ','').replace('"','')
				try:
					gcc=GC(seq)
					sys.stdout.write("%s\t%s\n"%(gene_id_bak,gcc))
				except NameError: 
					pass

				# initialize:
				seq=""
				gene_id_bak=gene_id

			if split_line[2]=='exon':
				chrom=split_line[0]
				start=int(split_line[3])-1
				end=int(split_line[4])
				seq+=genome[chrom][start:end].seq


else: 
	sys.stderr.write('3rd option must be either gene or exon')