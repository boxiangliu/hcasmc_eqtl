import matplotlib as mpl
mpl.use('Agg')
import pybedtools as bt
from os import listdir
import pandas as pd 
import matplotlib.pyplot as plt

IN_DIR='../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/'
OUT_DIR='../processed_data/hcasmc_specific_open_chromatin/intersect/'


in_files=listdir(IN_DIR)
bed_files=dict()
for i in in_files:
	epigenome=i.replace('.bed','')
	bed_files[epigenome]=bt.BedTool(IN_DIR+'/'+i)


op=bt.BedTool()

multi_inter=op.multi_intersect(i=[bed_files[x].fn for x in bed_files],cluster=True,header=True,names=[x for x in bed_files],output=OUT_DIR+'/intersect.bed')

x=pd.read_table(OUT_DIR+'/intersect.bed')
plt.hist(x.num,109)
plt.savefig('../figures/hcasmc_specific_open_chromatin/num_tissues_vs_peak_frequency.pdf')


