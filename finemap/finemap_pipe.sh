# Finemap: 
bash finemap/finemap/howson_preprocess.sh
Rscript finemap/finemap/howson_data.R
bash finemap/finemap/finemap.sh

Rscript finemap/finemap/rasqual_data.R
Rscript finemap/finemap/howson_rasqual_coloc.R

Rscript finemap/finemap/nikpay_rasqual_coloc.R

Rscript finemap/finemap/atacseq_overlap.R
Rscript finemap/finemap/atacseq_overlap.gviz.R

Rscript finemap/finemap/motif_match.R

Rscript finemap/finemap/tarid_tcf21_coloc.R

# Locuszoom:
Rscript finemap/locuszoom/locuscatter.R