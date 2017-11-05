import sys
for line in sys.stdin:
	line=line.strip()
	if not line.startswith('#'):
		split_line=line.split('\t')
		split_line[0]='chr'+split_line[0]
		split_line[3]='chr'+split_line[3]
		line='\t'.join(split_line)
	sys.stdout.write(line+'\n')