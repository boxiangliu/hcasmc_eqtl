#!/usr/bin/env python
# boxiang liu
# durga
# count the number of significant associations from matrix eqtl output
# examples: 
# python $scripts/count_num_sig_association.py /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/find_optimum_num_PEER_factors_matrixeqtl/pc3.peer1.2.cis.txt 0.1,0.05,0.01,0.001 fdr
# python $scripts/count_num_sig_association.py /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/find_optimum_num_PEER_factors_matrixeqtl/pc3.peer1.2.cis.txt 1e-3,1e-5,1e-8 pvalue

# modules:
import sys

# command line: 
input_file=sys.argv[1]
cutoff=sys.argv[2]
switch=sys.argv[3]

# output_file=sys.argv[3]


# convert cutoff into a list: 
cutoff=cutoff.strip().split(',')
cutoff=[float(x) for x in cutoff]
cutoff.sort()



# count the number of associations with FDR less than the cutoff: 
num_sig_association=[0 for x in cutoff]
col=4 if switch=='pvalue' else 5
with open(input_file,'r') as inp:
	for line in inp:
		split_line=line.strip().split()
		if split_line[0]=='SNP': continue # skip header
		stat=float(split_line[col])
		if stat>cutoff[-1]:
			break
		else:
			for i,c in enumerate(cutoff):
				is_sig=(stat<=c)
				num_sig_association[i]=num_sig_association[i]+int(is_sig)

# write output: 
num_sig_association=[str(x) for x in num_sig_association]
sys.stdout.write("\t".join(num_sig_association))
