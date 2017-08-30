import sys

fn='../data/encode/dnase_seq_2007_2012/controlled_vocabulary/cv.ra'

with open(fn,'r') as f:
	for line in f:
		if line.startswith('### Mouse Cell Lines'): break
		if line.startswith('#'): continue

		if line.strip()=='':
			if 'out' not in globals(): # header line
				out_line='term\ttag\ttype\torganism\tdescription\ttissue\tvendorName\tvendorId\torderUrl\tkaryotype\tlineage\ttermId\ttermUrl\tcolor\tsex\ttier\tprotocol\tcategory\n'
				sys.stdout.write(out_line)

			else: # data line
				out_line=''
				for field in ['term',
							  'tag',
							  'type',
							  'organism',
							  'description',
							  'tissue',
							  'vendorName',
							  'vendorId',
							  'orderUrl',
							  'karyotype',
							  'lineage',
							  'termId',
							  'termUrl',
							  'color',
							  'sex',
							  'tier',
							  'protocol',
							  'category']:
					try:
						out_line+=out[field]+'\t'
					except KeyError:
						out_line+='missing\t'

				out_line=out_line.strip()+'\n'
				sys.stdout.write(out_line)


			out=dict() # initialize or re-initialize 'out'

		else:
			field,annotation=line.strip().split(' ',1)
			out[field]=annotation

