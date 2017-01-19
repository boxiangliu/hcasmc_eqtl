# generate metasoft input:
#### modules: 
import sys
import scipy.stats as stats
import numpy as np
import gzip

#### constants: 
DEBUG=True
eqtl_file=sys.argv[1]
specificity_file=sys.argv[2]
NAfill=float(sys.argv[3])
# eqtl_file='../processed_data/eqtl_and_atacseq/eqtl.txt'
# specificity_file='../processed_data/eqtl_and_atacseq/specificity.txt'

TISSUE_ORDER=["Adipose_Subcutaneous",
"Adipose_Visceral_Omentum",
"Adrenal_Gland",
"Artery_Aorta",
"Artery_Coronary",
"Artery_Tibial",
"Brain_Anterior_cingulate_cortex_BA24",
"Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere",
"Brain_Cerebellum",
"Brain_Cortex",
"Brain_Frontal_Cortex_BA9",
"Brain_Hippocampus",
"Brain_Hypothalamus",
"Brain_Nucleus_accumbens_basal_ganglia",
"Brain_Putamen_basal_ganglia",
"Breast_Mammary_Tissue",
"Cells_EBV-transformed_lymphocytes",
"Cells_Transformed_fibroblasts",
"Colon_Sigmoid",
"Colon_Transverse",
"Esophagus_Gastroesophageal_Junction",
"Esophagus_Mucosa",
"Esophagus_Muscularis",
"HCASMC",
"Heart_Atrial_Appendage",
"Heart_Left_Ventricle",
"Liver",
"Lung",
"Muscle_Skeletal",
"Nerve_Tibial",
"Ovary",
"Pancreas",
"Pituitary",
"Prostate",
"Skin_Not_Sun_Exposed_Suprapubic",
"Skin_Sun_Exposed_Lower_leg",
"Small_Intestine_Terminal_Ileum",
"Spleen",
"Stomach",
"Testis",
"Thyroid",
"Uterus",
"Vagina",
"Whole_Blood"]

#### function:
def peek(x):
	if DEBUG:
		sys.stderr.write(str(x)+'\n')


def parse_line(line):
	split_line=line.strip().split()
	parsed_line={}
	parsed_line['ID']=split_line[0]
	parsed_line['pval']=float(split_line[1])
	parsed_line['beta']=float(split_line[2])
	parsed_line['se']=float(split_line[3])
	parsed_line['tissue']=split_line[4]
	return parsed_line


def format_output(association,tissue_order,ID):
	out_line=[ID]
	for tissue in tissue_order:
		if tissue in association:
			out_line=out_line+[str(association[tissue]['beta'])]
		else:
			out_line=out_line+["0"]
	joint_line="\t".join(out_line)
	return joint_line


def pass_pval_filter(association,p):
	passed=False
	pval=[value['pval'] for value in association.itervalues()]
	if min(pval)<p:
		passed=True
	return passed


def calc_specificity(beta,idx):
	beta=np.absolute(np.array(beta))
	p=beta/sum(beta)
	H=stats.entropy(p,base=2)
	Q=H-np.log2(p[idx])
	return Q

def calc_mean(association):
	beta=[v['beta'] for k,v in association.iteritems()]
	mean=float(sum(beta))/len(beta)
	return mean

def pass_filter(association):
	return ('HCASMC' in association) and (len(association) > 1)

def get_specificity_output(association,tissue_order,ID,idx,NAfill=-1):
	if NAfill == -1:
		NAfill=calc_mean(association)
	out_line=[ID]
	beta=[]
	for tissue in tissue_order:
		if tissue in association:
			beta+=[association[tissue]['beta']]
		else:
			beta+=[NAfill]
	Q=calc_specificity(beta, idx)
	out_line+=[str(Q)]
	return '\t'.join(out_line)


#### main 
association={}
output_line_ID=""
n_lines=0
n_input_lines=0
f_eqtl=gzip.open(eqtl_file,'w')
f_specificity=gzip.open(specificity_file,'w')

for line in sys.stdin:
	n_input_lines=n_input_lines+1
	parsed_line=parse_line(line)
	if parsed_line['ID']!=output_line_ID and output_line_ID!="":
		
		if pass_filter(association):
			# write current line to stdout:
			eqtl_output=format_output(association, TISSUE_ORDER, output_line_ID)
			f_eqtl.write(eqtl_output+'\n')

			# calculate specificity:
			specificity_output=get_specificity_output(association, TISSUE_ORDER, output_line_ID, TISSUE_ORDER.index('HCASMC'),NAfill=-1)
			f_specificity.write(specificity_output+'\n')

			# update number of written line: 
			n_lines=n_lines+1

		# progress report: 
		if n_lines%10000==0:
			sys.stderr.write(str(n_lines)+" output lines\n")

		# re-initialize association: 
		association={}

	# progress report: 
	if n_input_lines%10000000==0:
		sys.stderr.write(str(n_input_lines)+" input lines\n")


	# add to current line: 
	association[parsed_line['tissue']]={'pval':parsed_line['pval'],'beta':parsed_line['beta'],'se':parsed_line['se']}


	# update current ID:
	output_line_ID=parsed_line['ID']
