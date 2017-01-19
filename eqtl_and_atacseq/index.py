#!/usr/bin/env python 
import sys,pickle
fin=sys.argv[1]
field=int(sys.argv[2])
delim=sys.argv[3]

## Functions: 
def build_index(fin,delim='\t',field=1):
	index={}
	with open(fin, 'r') as f:
		eof=False
		byte=0
		while not eof:
			line=f.readline()
			eof=(line=="")
			split_line=line.strip().split(delim)
			try:
				key=split_line[field-1]
			except IndexError:
				continue
			if key not in index:
				index[key]=[byte]
			else:
				index[key].append(byte)
			byte=f.tell()
	return index

def query(fin,index,key):
	with open(fin,'r') as f:
		for pos in index[key]:
			f.seek(pos)
			sys.stdout.write(f.readline())


## Main:
# Build index: 
index=build_index(fin,delim=delim,field=field)
pickle.dump(index, open(fin+'.idx','w'))


# Query: 
fin='iris.data'
index=pickle.load(open('iris.data.idx','rb'))
key='5.2'
query(fin, index, key)



