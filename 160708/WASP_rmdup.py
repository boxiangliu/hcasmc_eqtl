import sys 
sys.path.append('/srv/persistent/bliu2/HCASMC_eQTL/scripts/')
from utils import *
from multiprocessing import Pool

def WASPrmdup(bam, paired_end = True):
	WASP_mapping = '/srv/persistent/bliu2/tools/WASP-0.1/mapping'


	output = bam.replace('sorted.bam','rmdup.bam')
	if paired_end: 
		cmd = 'python %s/rmdup_pe.py %s %s'%(WASP_mapping, bam, output)
	else: 
		cmd = 'python %s/rmdup.py %s %s'%(WASP_mapping, bam, output)
	call(cmd)
	sort(output)
	index(output.replace('bam', 'sorted.bam'))

	print "[ WASPrmdup ] %s finished."%bam


if __name__ == '__main__':
	setwd(sys.argv[1])

	sample_list = sys.argv[2]
	sample_dirs = readSampleList(sample_list)

	keep_bams = ["%s/%s"%(sample_dir,'wasp.keep.merged.sorted.bam') for sample_dir in sample_dirs]
	
	pool = Pool(processes = 10)
	pool.map(WASPrmdup, keep_bams)
	
	reportDone()