# 16/06/03:
# setup: 
# prepare GTEx and hcasmc rpkm files (all genes): 
mkdir $processed_data/rpkm
bash $scripts/prepare_gtex_rpkm.sh $processed_data/rpkm
ln ../processed_data/160527/combined.rpkm $processed_data/rpkm/hcasmc.rpkm


# hierarchical clustering:
Rscript $scripts/hclust.R \
	-rpkm=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm \
	-coldata=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col \
	-figure=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/hclust.pdf

# multidimensional scaling (2D):
Rscript $scripts/mds.R \
	-rpkm=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm \
	-coldata=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col \
	-tissue_names=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.txt \
	-figure=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/

# multidimensional scaling (3D).
Since rgl is not installed on durga, the script needs to be run locally. 
Rscript $scripts/mds.3D.R \
	-rpkm=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.rpkm \
	-coldata=/srv/persistent/bliu2/HCASMC_eQTL/processed_data/160603/combined.col \
	-tissue_names=/srv/persistent/bliu2/HCASMC_eQTL/scripts/160603/collapsed_tissue_names.txt \
	-figure=/srv/persistent/bliu2/HCASMC_eQTL/figures/160603/



