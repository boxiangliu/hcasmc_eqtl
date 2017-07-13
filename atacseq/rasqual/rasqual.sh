# run RASQUAL
# bash rasqual.sh <input parameter file> <line number> Y.bin K.bin X.bin VCF

param_file=$1
line_num=$2
Y=$3
K=$4
vcf_file=$5
out_dir=$6


param=($(cat $param_file | sed "${line_num}q;d"))
gene_id=${param[0]}
gene_name=${param[1]}
region=${param[2]}
n_rsnp=${param[3]}
n_fsnp=${param[4]}
exon_start_positions=${param[5]}
exon_end_positions=${param[6]}
feat_id=$(grep $gene_id -n ../processed_data/atacseq/rasqual//expression/Y.txt | cut -d":" -f1,1)
window_size=200000
n_sample=8
echo id: $gene_id 
echo name: $gene_name 
echo region: $region
echo reference snps: $n_rsnp
echo feature snps: $n_fsnp
echo feature id: $feat_id
echo Y: $Y
echo K: $K
echo vcf: $vcf_file

if [[ -s $out_dir/${gene_id}.txt ]]; then 
	echo $out_dir/${gene_id}.txt exist! skipping...
else  
	tabix $vcf_file $region | \
	/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
	-y $Y -k $K \
	-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
    -s $exon_start_positions -e $exon_end_positions \
    --imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
    --cis-window-size $window_size \
    -f $gene_id --n_threads 1 \
    --force > $out_dir/${gene_id}.txt
fi 
