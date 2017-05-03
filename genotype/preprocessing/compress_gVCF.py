from utils import *


setwd(sys.argv[1])
sample_dirs = readSampleList(sys.argv[2])
pool = Pool(processes=len(sample_dirs))
pool.map(compressVcf, sample_dirs)
pool.map(indexVcf, sample_dirs)

reportDone()
