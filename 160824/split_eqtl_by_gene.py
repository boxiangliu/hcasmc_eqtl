# split the eqtl by gene id:
import sys
out_dir=sys.argv[1]
pid=""
out_string=""
n=0
for line in sys.stdin:
	split_line=line.strip().split()
	cid=split_line[0]
	# when end of current gene is reached:
	if cid!=pid and pid!="":
		out_file=out_dir+"/"+pid+'.txt'
		with open(out_file,'w') as out:
			out.write(out_string)
		out_string=""
		sys.stderr.write(pid+'\n')
	pid=cid	
	out_string=out_string+line

# when EOF is reached:
out_file=out_dir+"/"+cid+'.txt'
with open(out_file,'w') as out:
	out.write(out_string)
