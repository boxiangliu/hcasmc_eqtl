#!/bin/bash
# Nathan Abell
# Script to execute GSEA2 Pre-Ranked on a RNK-format file
# Usage: bash GSEAPreRanked.sh <RNK_file> <output_directory>


# Variable declarations:
export GSEA="/users/bliu2/tools/gsea/"

# Fucntions:
GSEAPreranked(){

RNK=$1
OUTDIR=$2
LABEL=$3
LOG=$OUTDIR/GSEA_run.log

java -Xmx4G -cp $GSEA/gsea-3.0.jar \
xtools.gsea.GseaPreranked \
-rnd_seed 2017 \
-rnk $RNK \
-gmx $GSEA/MSigDB/h.all.v6.1.symbols.gmt \
-out $OUTDIR \
-rpt_label  $LABEL \
-plot_top_x 50 \
-collapse false > $LOG 2>&1

}

export -f GSEAPreranked


# Main: 
parallel -j3 GSEAPreranked \
../processed_data/differential_expression/DESeq2/{}.rnk \
../processed_data/differential_expression/GSEA/ \
{} ::: meta_artery meta_heart meta_fibroblast


# Make table Name (NES,FWER p-value): 
# awk '{FS="\t";OFS=""}{printf "%s (%.02f,%.02f)\n",$1,$6,$8}' ../processed_data/differential_expression/GSEA/meta_fibroblast.GseaPreranked.*/gsea_report_for_na_pos_*.xls
# awk '{FS="\t";OFS=""}{printf "%s (%.02f,%.02f)\n",$1,$6,$8}' ../processed_data/differential_expression/GSEA/meta_artery.GseaPreranked.*/gsea_report_for_na_pos_*.xls
# awk '{FS="\t";OFS=""}{printf "%s (%.02f,%.02f)\n",$1,$6,$8}' ../processed_data/differential_expression/GSEA/meta_heart.GseaPreranked.*/gsea_report_for_na_pos_*.xls
