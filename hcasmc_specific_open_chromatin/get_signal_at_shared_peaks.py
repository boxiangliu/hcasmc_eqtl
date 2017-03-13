import pybedtools as bt
from os import listdir
from os import remove
import pandas as pd 
import atexit


# Read data: 
x=pd.read_table('../processed_data/hcasmc_specific_open_chromatin/intersect/intersect.bed')


# Creating id column:
x['id']=x[['chrom','start','end','num']].apply(lambda x: '_'.join([str(y) for y in x]),axis=1)


# Rearrange columns: 
cols=x.columns.tolist()
cols=cols[0:5]+[cols[len(cols)-1]]+cols[6:(len(cols)-2)]
x=x[cols]


# Rename columns: 
x.rename(columns={'start':'chromStart','end':'chromEnd'},inplace=True)


# Output to temporary file: 
tmp_file='../processed_data/hcasmc_specific_open_chromatin/intersect/tmp_file.bed'
x.iloc[:,0:6].to_csv(tmp_file,sep='\t',header=False,index=False)
atexit.register(remove,tmp_file)


# Create bedtool object:
intersect=bt.BedTool(tmp_file)


# Get signal at shared peaks: 
IN_DIR='../processed_data/hcasmc_specific_open_chromatin/encode_plus_hcasmc_filt/'
in_files=listdir(IN_DIR)
op=bt.BedTool()


for i in in_files:
	print(i)
	a=bt.BedTool(IN_DIR+'/'+i)
	res=op.intersect(a=a,b=intersect,wa=True,wb=True)


	id=list()
	signalValue=list()
	for interval in res:
		id.append(interval[9])
		signalValue.append(interval[3])


	sample=i.replace('.bed','')
	df=pd.DataFrame({'id':id,sample:signalValue})


	df=df[['id',sample]]
	df.to_csv('../processed_data/hcasmc_specific_open_chromatin/intersect/samplewise_count/intersect.%s.count'%sample,sep='\t',header=True,index=False)
