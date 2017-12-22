# Subsample GTEx to 52 and map eQTL
# Boxiang Liu
# 2017-12-21

# subset to covariates to first 10: 
bash $scripts/160816/subset_covariates.sh 


# subsample tissues: 
Rscript $scripts/160816/subsample.lists.by.tissue.R


# run fastqtl nominal pass on subsampled GTEx tissues: 
bash $scripts/160816/nominal.pass.subsamples.v6p.sh


# run fastqtl permutation pass on subsampled GTEx tissues: 
bash $scripts/160816/permutation.pass.subsamples.v6p.sh


#### replication of HCASMC eQTLs in GTEx tissues: 
# select HCASMC eQTLs with FDR < 0.05: 
cat ../processed_data/160805/hcasmc.eqtl.pc4.peer8.perm.padj.txt | \
awk 'BEGIN{OFS="\t"} {if ($19<=0.05) {print $1,$6"_b37"}}' | \
sed "s/chr//" > ../processed_data/160816/hcasmc.eqtl.pc4.peer8.perm.padj.fdr05.txt


# subset GTEx tissue to association significant in HCASMC:
mkdir $processed_data/160816/replication/
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
python $scripts/160816/subset_eQTLs.py \
../processed_data/160816/hcasmc.eqtl.pc4.peer8.perm.padj.fdr05.txt \
$data/gtex/v6p/v6p_fastQTL_allpairs_FOR_QC_ONLY/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz \
../processed_data/160816/replication/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.sig_in_HCASMC.txt


# calculate pi1 statistics for GTEx tissues
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
Rscript $scripts/160816/pi1_calc.R \
../processed_data/160816/replication/{}_Analysis.v6p.FOR_QC_ONLY.allpairs.sig_in_HCASMC.txt \
{} > ../processed_data/160816/replication/pi1.txt


# plot pi1 statistic:
Rscript $scripts/160816/pi1_plot.R \
../processed_data/160816/replication/pi1.txt \
../figures/160816/pi1.pdf


# subset subsampled GTEx tissue to association significant in HCASMC:
mkdir $processed_data/160816/replication_subsample/
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
python $scripts/160816/subset_eQTLs.py \
../processed_data/160816/hcasmc.eqtl.pc4.peer8.perm.padj.fdr05.txt \
$processed_data/160816/subsampling/{}/{}_52.allpairs.txt.gz \
../processed_data/160816/replication_subsample/{}_52.allpairs.sig_in_HCASMC.txt


# calculate pi1 statistics for GTEx tissues
parallel -a $data/gtex/gtex.v6p.eqtl.tissues.txt -j44 \
Rscript $scripts/160816/pi1_calc.R \
../processed_data/160816/replication_subsample/{}_52.allpairs.sig_in_HCASMC.txt \
{} > ../processed_data/160816/replication_subsample/pi1.txt


# plot pi1 statistic:
Rscript $scripts/160816/pi1_plot.R \
../processed_data/160816/replication_subsample/pi1.txt \
../figures/160816/pi1_subsample.pdf

