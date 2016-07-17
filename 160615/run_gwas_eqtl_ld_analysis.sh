#!/bin/bash 

# command line args: 
fdr=$1
eqtldir=$2
alpha=$3
lddir=$4
# fdr=0.05
# alpha=1e-8
# eqtldir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/compare_hcasmc_and_gtex_fdr${fdr}
# lddir=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/LD_fdr${fdr}_alpha${alpha}


# constants: 
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160615
processed_data=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615
figures=/srv/persistent/bliu2/HCASMC_eQTL/figures/160615



# mkdir:
if [[ ! -d $eqtldir ]]; then mkdir $eqtldir; fi 
if [[ ! -d $lddir ]]; then mkdir $lddir; fi 


# get top eQTLs from HCASMC:
echo "getting eqtls..."
Rscript $scripts/get_top_eqtls.hcasmc.R \
-input=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160530/find_optimum_num_PEER_factors/fastqtl.pc3.peer8.padj.txt \
-output=$eqtldir/HCASMC.txt \
-fdr=$fdr


# get the number of significant HCASMC eQTLs: 
num_sig=$(wc -l $eqtldir/HCASMC.txt | awk '{print $1}')


# get top eQTLs from GTEx tissue:
bash $scripts/run_get_top_eqtls.gtex.sh $num_sig $eqtldir


# get top GWAS hits:
echo "getting gwas hits..."
Rscript $scripts/get_top_gwas.R \
-input=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160615/gwas/cad.add.160614.website.txt \
-output=$eqtldir/GWAS.txt \
-alpha=$alpha


# calculate the LD between GWAS and eQTL hits:
echo "calculating ld..."
bash $scripts/calculate_ld.sh $eqtldir $lddir


# plot the number of significant LD variants:
echo "plotting..."
cat $eqtldir/*.txt > $lddir/combined.txt
Rscript $scripts/plot_ld_by_tissue.R \
	-ld=$lddir/target_variants.ld \
	-bim=$lddir/target_variants.bim \
	-variants=$lddir/combined.txt \
	-figure=$figures/ld_distribution_by_tissue.fdr$fdr.alpha$alpha.pdf