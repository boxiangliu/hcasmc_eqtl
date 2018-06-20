#  -*- sh -*-
#$ -S /bin/bash
#$ -q res.q
#$ -cwd
#$ -V
#$ -l mem=12G
#$ -pe smp 4
##$ -m e
##$ -M skim@nygenome.org
#$ -j y
#$ -N  xcell
set -x

## usage: qsub -p -10 ./Run_xCell_custom.sig.pancreas_vm.sh 

module load R/3.2.2

dat=v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY_xCell.gct.gz # needs about 10min to read, and 6h to compute scores
name=xCell_rpkm_res_all.tissues_pancreas 
core=4 
signature=pancreas_all.markers.Rda


Rscript Run_xCell_custom.sig.pancreas_vm.R ${dat} ${name} ${core} ${signature}
