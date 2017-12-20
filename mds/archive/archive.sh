# prepare GTEx and hcasmc rpkm files (all genes): 
mkdir $processed_data/rpkm
bash $scripts/prepare_gtex_rpkm.sh $processed_data/rpkm
ln ../processed_data/160527/combined.rpkm $processed_data/rpkm/hcasmc.rpkm

