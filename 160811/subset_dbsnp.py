# subset to lines in dbSNP based on input file s:
import sys

# constants: 
DBSNP='/srv/persistent/bliu2/shared/dbsnp/snp146.txt'


# read file; will subset DBSNP to this file: 
s=set(["_".join(line.strip().split()) for line in sys.stdin])
sys.stderr.write(str(len(s))+' items\n')


# subset dbSNP to set s:
n_in_lines=0
n_out_lines=0
with open(DBSNP,'r') as f:
	for line in f:
		n_in_lines=n_in_lines+1
		split_line=line.strip().split()
		chr_pos="_".join([split_line[1],split_line[3]])
 
		# check whether line is in set s: 
		if chr_pos in s:
			n_out_lines=n_out_lines+1
			out_line="\t".join([split_line[1],split_line[3],split_line[4],split_line[6],split_line[9]])
			print out_line

			# progress report: 
			if n_out_lines%100==0:
				sys.stderr.write(str(n_out_lines)+" output lines\n")

		# progress report: 
		if n_in_lines%1000000==0:
			sys.stderr.write(str(n_in_lines)+" input lines\n")