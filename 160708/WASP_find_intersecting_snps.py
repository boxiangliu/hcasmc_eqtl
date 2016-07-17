# run find_intersecting_snps.py: 
# TODO: save the output of find_intersecting_snps.py to a log file.
 
import sys
sys.path.append('/srv/persistent/bliu2/HCASMC_eQTL/scripts')
from utils import *
from multiprocessing import Pool


def copy(directory, src='Aligned.out.sorted.rg.uniq.bam', dst='wasp.bam'):
	''' copy directory/src to directory/dst'''
	src = directory + '/' + src
	dst = directory + '/' + dst
	if not os.path.exists(dst): 
		print "[ copy ] copying %s -> %s"%(src, dst)
		shutil.copyfile(src, dst)
	return dst

def WASPfindIntersectingSNPs(bam, SNP_file_directory='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160708/WASP_SNPs/', paired_end = True):
	''' param str bam: path to bam file
		param str SNP_file_directory: path to SNP file directory

		return: none 
		note: Output will be in the same directory as bams'''
	WASP_mapping = '/srv/persistent/bliu2/tools/WASP/mapping'
	if paired_end: 
		mode = '-p'
	else: 
		mode = ''

	cmd = 'python %s/find_intersecting_snps.py %s %s %s'%(WASP_mapping, mode, bam, SNP_file_directory)
	call(cmd)

	print "[ WASPfindIntersectingSNPs ] %s finished."%bam

if __name__ == '__main__':
	setwd(sys.argv[1])

	sample_list = sys.argv[2]
	sample_dirs = readSampleList(sample_list)

	# find intersecting snps: 
	pool = Pool(processes = len(sample_dirs))
	bams = pool.map(copy, sample_dirs, 4) # from secondPass.Aligned.sortedByCoord.out.bam to wasp.bam
	pool.map(WASPfindIntersectingSNPs, bams, 4)

	reportDone()