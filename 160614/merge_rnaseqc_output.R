#!/usr/bin/env Rscript
# boxiang liu
# durga
# merge RNA-seQC output

args=commandArgs(T,T)
input_list_file=args$input
output_file=args$output
# input_list_file='160614/merge_rnaseqc_output.eqtl.txt'
# output_file='/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160614/rnaseqc.hcasmc_eqtl.reads.gct'

input_list=scan(input_list_file,what=character())
output=data.table()

for (input_file in input_list){
	input=fread(input_file,skip=2,header=T)
	if (ncol(output)==0){
		output=input
	} else {
		stopifnot(output$Name==input$Name)
		output=data.table(output,input[,3,with=F])
	}
}

writeLines(sprintf('#1.2\n%s\t%s',nrow(output),ncol(output)-2),output_file)
write.table(output,output_file,append=T,quote=F,sep='\t',row.names=F,col.names=T)