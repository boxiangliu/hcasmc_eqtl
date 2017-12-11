gene_id=$1
eqtl_fn=$2
out_fn=$3

echo INFO - gene: $gene_id
echo INFO - eqtl: $eqtl_fn
echo INFO - output: $out_fn

out_dir=$(dirname $out_fn)
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

zcat $eqtl_fn | grep $gene_id > $out_fn.tmp

Rscript finemap/locuszoom/find_rsid.R $out_fn.tmp $out_fn
rm $out_fn.tmp
