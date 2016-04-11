from utils import *

setwd(sys.argv[1])
sample_dirs = readSampleList(sys.argv[2])

for sample_dir in sample_dirs:
	sample_dir = sample_dir.strip()
	cmd = "tabix %s/raw_variants.g.vcf.gz chrY:59373566"%sample_dir
	call(cmd)

reportDone()