#### modules: 
import sys


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

#### function:
def parse_line(line):
	split_line=line.strip().split()
	parsed_line={}
	parsed_line['ID']=split_line[0]
	parsed_line['pval']=float(split_line[1])
	parsed_line['beta']=float(split_line[2])
	parsed_line['se']=float(split_line[3])
	parsed_line['tissue']=split_line[4]
	return parsed_line

def id_changed(this_id, last_id):
	changed=((this_id!=last_id) and (last_id!=""))
	return changed

def format_output(association,tissue_order,ID):
	out_line=''
	for tissue in TISSUE_ORDER:
		if tissue in association:
			out_line=out_line+ID+'\t'+'\t'.join([str(association[tissue][k]) for k in ['pval','beta','se']])+'\t'+tissue+'\n'
	return out_line


def pass_filter(association,num_tissue=10,sig_pval=1e-6,noise_pval=5e-2):
	if ('HCASMC' in association) and (len(association) >= num_tissue): 
		is_noise=[association[tissue]['pval']>noise_pval for tissue in association if association[tissue] != 'HCASMC']
		all_noise=reduce(lambda x,y: x and y, is_noise)
		return (association['HCASMC']['pval']<sig_pval) and (all_noise)
	else: 
		return False




#### main 
association={}
output_line_ID=""
n_lines=0
n_input_lines=0


for line in sys.stdin:
	n_input_lines=n_input_lines+1
	parsed_line=parse_line(line)
	if id_changed(parsed_line['ID'],output_line_ID):
		if pass_filter(association,num_tissue=10,sig_pval=1e-6,noise_pval=5e-2):
			# write current line to stdout:
			out_line=format_output(association, TISSUE_ORDER, output_line_ID)
			sys.stdout.write(out_line)

			# update number of written line: 
			n_lines=n_lines+1

		# re-initialize association: 
		association={}

	# progress report: 
	if n_input_lines%10000000==0:
		sys.stderr.write(str(n_input_lines)+" input lines\n")


	# add to current line: 
	association[parsed_line['tissue']]={'pval':parsed_line['pval'],'beta':parsed_line['beta'],'se':parsed_line['se']}


	# update current ID:
	output_line_ID=parsed_line['ID']

# output the last line:
if pass_filter(association):
	# write current line to stdout:
	out_line=format_output(association, TISSUE_ORDER, output_line_ID)
	sys.stdout.write(out_line)

