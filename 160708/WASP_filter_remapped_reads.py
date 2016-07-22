from multiprocessing import Pool
import sys
sys.path.append('/srv/persistent/bliu2/HCASMC_eQTL/scripts')
from utils import *

def WASPfilterRemappedReads(params, paired_end = True):
	(to_remap_bam, remapped_bam) = params 

	WASP_mapping = '/srv/persistent/bliu2/tools/WASP/mapping'

	if paired_end: 
		mode = '-p'
	else: 
		mode = ''

	output_bam = to_remap_bam.replace('to.remap.bam','remap.keep.bam')
	to_remap_num = to_remap_bam.replace('to.remap.bam',"to.remap.num.gz")
	cmd = "python %s/filter_remapped_reads.py %s %s %s %s %s"%(WASP_mapping, mode, to_remap_bam, remapped_bam, output_bam, to_remap_num)
	call(cmd)

	# merge, sort and index:		
	keep_bam = to_remap_bam.replace('to.remap.bam','keep.bam')
	remap_keep_bam = to_remap_bam.replace('to.remap.bam','remap.keep.bam')
	merged_bam =  to_remap_bam.replace('to.remap.bam','keep.merged.bam')
	merge([keep_bam, remap_keep_bam], merged_bam)
	
	sort(merged_bam)

	sorted_bam = merged_bam.replace('bam','sorted.bam')
	index(sorted_bam)

	print "[ WASPfilterRemappedReads ] %s finished."%(os.path.basename(to_remap_bam))



if __name__ == '__main__':
	setwd(sys.argv[1])

	sample_list = sys.argv[2]
	sample_dirs = readSampleList(sample_list)

	# filter remapped reads:
	params = [("%s/%s"%(sample_dir,'wasp.to.remap.bam'),"%s/%s"%(sample_dir,'wasp.remapped.Aligned.out.bam')) for sample_dir in sample_dirs]
	pool = Pool(processes = len(sample_dirs))
	pool.map(WASPfilterRemappedReads, params, 3)

	reportDone()