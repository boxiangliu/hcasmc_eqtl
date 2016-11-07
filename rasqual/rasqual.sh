# run RASQUAL
# bash rasqual.sh <input parameter file> <line number> Y.bin K.bin X.bin VCF

param_file=$1
line_num=$2
Y=$3
K=$4
X=$5
vcf_file=$6

param=($(cat $1 | sed "${line_num}q;d"))
gene_id=${param[0]}
gene_name=${param[1]}
region=${param[2]}
n_rsnp=${param[3]}
n_fsnp=${param[4]}
exon_start_positions=${param[5]}
exon_end_positions=${param[6]}
feat_id=$(grep $gene_id -n /srv/persistent/bliu2/HCASMC_eQTL/processed_data//rasqual/Y.tidy.txt | cut -d":" -f1,1)
window_size=2000000
n_sample=52
echo id: $gene_id 
echo name: $gene_name 
echo region: $region
echo reference snps: $n_rsnp
echo feature snps: $n_fsnp
echo feature id: $feat_id

if [[ -e ../processed_data/rasqual/output/${gene_id}_${gene_name}.txt ]]; then 
	echo "../processed_data/rasqual/output/${gene_id}_${gene_name}.txt exist! skipping..."
else  
	tabix $vcf_file $region | \
	/srv/persistent/bliu2/tools/rasqual/src/rasqual \
	-y $Y \
	-k $K \
	-x $X \
	-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
    -s $exon_start_positions -e $exon_end_positions \
    --imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
    -f $gene_name --n_threads 1 \
    --force -v > ../processed_data/rasqual/output/${gene_id}_${gene_name}.txt
fi 
