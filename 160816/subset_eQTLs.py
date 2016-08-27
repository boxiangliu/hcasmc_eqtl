import sys,gzip

# command line argument: 
filter_file=sys.argv[1]
eqtl_file=sys.argv[2]
out_file=sys.argv[3]


# read filter file:
filt=set()
n_filt_line=0
with open(filter_file,'r') as f:
	for line in f:
		n_filt_line=n_filt_line+1
		split_line=line.strip().split()
		filt.add((split_line[0],split_line[1]))

assert len(filt)==n_filt_line, 'entries in %s file are not unique'%filter_file


# subset eQTL file to those in the filter file:
n_output_line=0
with gzip.open(eqtl_file,'r') as f, open(out_file,'w') as o:
	for line in f: 
		split_line=line.strip().split()
		pheno=split_line[0]
		geno=split_line[1]
		if (pheno,geno) in filt:
			n_output_line=n_output_line+1
			o.write(line)

sys.stderr.write('%s lines written\n'%str(n_output_line))

