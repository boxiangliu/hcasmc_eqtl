import sys

min_pval=2
fid=""
sid=""
for line in sys.stdin:
	split_line=line.strip().split()

	if fid=="" and sid=="":
		fid,sid=split_line[0:2]

	if split_line[0]==fid and split_line[1]==sid:

		if float(split_line[3]) < min_pval:
			kept_line=line

	else:
		sys.stdout.write(kept_line)

		min_pval=float(split_line[3])
		kept_line=line

	fid,sid=split_line[0],split_line[1]

sys.stdout.write(kept_line)
