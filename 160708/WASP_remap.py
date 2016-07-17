from utils import *
from multiprocessing import Pool
from functools import partial
from align import removeGenome, loadGenome
debug = False
def WASPremap(read_pair, splice_junctions, genome_dir, n_thread = 8):
	'''Output will be in the same directory as reads'''

	read1 = read_pair[0]
	read2 = read_pair[1]
	prefix = read1.replace('remap.fq1.gz','remapped.')

	cmd = "STAR runThreadN %s --genomeDir %s --genomeLoad LoadAndKeep --readFilesIn %s %s --readFilesCommand zcat --sjdbFileChrStartEnd %s --outSAMtype BAM Unsorted --outFileNamePrefix %s"%(str(n_thread), genome_dir, read1, read2, splice_junctions, prefix)
	call(cmd)
		
	print "[ WASPremap ] %s finished."%os.path.dirname(read_pair[0])

if __name__ == '__main__':
	setwd(sys.argv[1])

	sample_list = sys.argv[2]
	sample_dirs = readSampleList(sample_list)

	# remap: 
	read_pairs = [("%s/%s"%(sample_dir,'wasp.remap.fq1.gz'), "%s/%s"%(sample_dir,'wasp.remap.fq2.gz')) for sample_dir in sample_dirs]
	if debug: print "read pairs"; print read_pairs
	
	splice_junctions = ' '.join(["%s/%s"%(sample_dir,'firstPass.SJ.out.tab') for sample_dir in ["ERR440421","ERR440422","ERR440423","ERR440424","ERR440425","ERR440426","ERR440427","ERR440428","ERR440429","ERR440430","ERR440431","ERR440432","ERR440433","ERR440434"]])
	# splice_junctions = ' '.join(["%s/%s"%(sample_dir,'firstPass.SJ.out.tab') for sample_dir in sample_dirs])
	if debug: print "splice junctions"; print splice_junctions
	
	genome_dir = '../../GRCm38/STAR_index_overhang75/'
	
	removeGenome(genome_dir)
	loadGenome(genome_dir)
	pool = Pool(processes = len(sample_dirs))
	func = partial(WASPremap, splice_junctions = splice_junctions, genome_dir = genome_dir)
	pool.map(func, read_pairs)
	removeGenome(genome_dir)

	reportDone()