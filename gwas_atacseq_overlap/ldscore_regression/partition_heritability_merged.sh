unset DISPLAY XAUTHORITY

export sumstats_fn='../processed_data/gwas_gene_overlap/ldscore_regression/gwas_sumstats/cad.sumstats.gz'
export baseline_annotation_dir='/srv/persistent/bliu2/shared/ldscore/baseline/'
export weight_dir='/srv/persistent/bliu2/shared/ldscore/weights_hm3_no_hla/'
export frq_dir='/srv/persistent/bliu2/shared/ldscore/1000G_frq/'

tissue_specific_annotation_dir=${1-'../processed_data/gwas_atacseq_overlap/ldscore_regression/ldscore_merged/'}
out_dir=${2-'../processed_data/gwas_atacseq_overlap/ldscore_regression/partition_heritability_merged/'}

mkdir -p $out_dir

partition_heritability(){
mkdir -p $2
echo INFO - in_dir: $1
echo INFO - out_dir: $2
if [[ -e "$2/cad.merged.results" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi 

python ~/tools/ldsc/ldsc.py \
--h2 $sumstats_fn \
--ref-ld-chr $1/merged.,$baseline_annotation_dir/baseline. \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $2/cad.merged \
--print-coefficients
}

export -f partition_heritability

partition_heritability_nobaseline(){

mkdir -p $2
echo INFO - in_dir: $1
echo INFO - out_dir: $2

if [[ -e "$2/cad.merged.nobaseline.results" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi 

python ~/tools/ldsc/ldsc.py \
--h2 $sumstats_fn \
--ref-ld-chr $1/merged. \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $2/cad.merged.nobaseline \
--print-coefficients
}

export -f partition_heritability_nobaseline


# for dir1 in jaccard_similarity_0.3 jaccard_similarity_0.4 jaccard_similarity_0.5 all_tissue; do
for dir1 in jaccard_similarity_0.3; do

echo INFO - $dir1
partition_heritability $tissue_specific_annotation_dir/$dir1/ $out_dir/$dir1/
partition_heritability_nobaseline $tissue_specific_annotation_dir/$dir1/ $out_dir/$dir1/

done