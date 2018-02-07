# generate metasoft input:

#### modules: 
import sys

#### constants: 
DEBUG=True

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
			out_line=out_line+[str(association[tissue]['beta']),str(association[tissue]['se'])]
		else:
			out_line=out_line+["NA","NA"]
	joint_line="\t".join(out_line)
	return joint_line


def pass_pval_filter(association,p):
	passed=False
	pval=[value['pval'] for value in association.itervalues()]
	if min(pval)<p:
		passed=True
	return passed


#### constants: 
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

#### main 
association={}
output_line_ID=""
n_lines=0
n_input_lines=0
if len(sys.argv) > 1:
	p=float(sys.argv[1])
else: 
	p=1e-3
sys.stderr.write('threshold: '+str(p)+'\n')
for line in sys.stdin:
	n_input_lines=n_input_lines+1
	parsed_line=parse_line(line)
	if parsed_line['ID']!=output_line_ID and output_line_ID!="":
		
		# write current line to stdout:
		passed=pass_pval_filter(association, p)
		if passed:
			output=format_output(association, TISSUE_ORDER, output_line_ID)
			print(output)


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
