import sys
# template="wgs_pipeline.scg3.sh"
template=sys.argv[1]
print template
samples="/mnt/data/WGS_HCASMC/WGS_HCASMC_fastq_list"
from pdb import set_trace 
with open(samples, 'r') as f:
	for line in f.readlines():

		line = line.strip()
		# print line
		if line.startswith("#"): continue
		if line == "": 
			sample = ""
			reads1 = ""
			reads2 = ""
			continue
		if ":" in line:
			sample = line.rstrip(":").lstrip("./")
		if sample in line or sample == "2109" or sample == "CA1401": 
			if "R1" in line: 
				reads1 = line
			elif "R2" in line: 
				reads2 = line
			else:
				pass 
		if len(sample) > 0 and len(reads1) > 0 and len(reads2) > 0: 
			output = template.replace('sh', '%s.sh'%sample)
			with open(template, 'r') as t, open(output, 'w') as output:
				for tline in t.readlines():
					if "placeholder1" in tline:
						tline = tline.replace("placeholder1", sample)
					elif "placeholder2" in tline: 
						tline = tline.replace('placeholder2', reads1)
					elif "placeholder3" in tline: 
						tline = tline.replace('placeholder3', reads2)
					else: 
						pass
					output.write(tline)
