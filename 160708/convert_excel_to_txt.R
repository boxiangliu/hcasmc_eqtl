# convert excel file to a txt file 
# subsettting a specific worksheet and columns

# library
library('XLConnect')


# command args:
args=commandArgs(T,T)
excel_file=args$excel
sheet=as.numeric(args$sheet)
subset=strsplit(args$subset,',')[[1]]
print(subset)
output_file=args$output

# excel_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rna_wgs_match.reduced_050616.xlsx'
# sheet=5
# subset=c('Path_to_alignment_on_valk','Path_to_alignment_on_durga')
# output_file='WASP_remap.sample_list.4.txt'

excel=readWorksheet(loadWorkbook(excel_file),sheet=sheet)
output=excel[,subset]

write.table(output,file=output_file,quote=F,row.names=F,sep='\t')